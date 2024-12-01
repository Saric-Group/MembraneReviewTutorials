from collections.abc import Callable
from dataclasses import dataclass
from math import nan
from pathlib import Path
from typing import Protocol, cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize
import scipy.special
import xarray as xr
from equilibration import series2ufloat
from numpy.typing import ArrayLike
from ovito.data import DataCollection, DataTable
from ovito.io import import_file
from ovito.modifiers import (
    ClearSelectionModifier,
    ClusterAnalysisModifier,
    ExpressionSelectionModifier,
    LoadTrajectoryModifier,
)
from tqdm.auto import tqdm
from z2u import zL2u


def modify_clusters(frame: int, data: DataCollection, threshold_count: int = 50):
    cluster_id, cluster_count = data.tables["clusters"].xy().T

    particle_cluster = data.particles_["Cluster_"]

    big_cluster_count = (cluster_count > threshold_count).sum()

    particle_cluster[:] = np.where(
        particle_cluster <= big_cluster_count, particle_cluster, 0
    )
    data.tables["clusters_"].delete_elements(cluster_id > big_cluster_count)

    particle_mol = np.array(data.particles["Molecule Identifier"])
    # don't assume all mols are prent
    mol_u, particle_mol_i = np.unique(particle_mol, return_inverse=True)

    particle_cluster = np.array(data.particles["Cluster"])
    mol_u_cluster = np.zeros(len(mol_u))
    particle_type = np.array(data.particles.particle_types)

    def min_cluster_per_mol(mol, cluster, m, c):
        result = np.full(m, c + 1)
        np.minimum.at(result, mol, cluster)
        return result

    sel = particle_type == 2
    v = min_cluster_per_mol(
        particle_mol_i[sel],
        particle_cluster[sel],
        len(mol_u),
        big_cluster_count + 1,
    )
    sel = v <= big_cluster_count
    mol_u_cluster[sel] = v[sel]

    particle_cluster = mol_u_cluster[particle_mol_i]

    data.particles_["Cluster_"][:] = particle_cluster

    cluster_id, cluster_count = np.unique(particle_cluster, return_counts=True)
    del data.tables["clusters"]

    key = "clusters"
    t = DataTable(title=key)
    t.x = t.create_property("id", data=cluster_id)
    t.y = t.create_property("count", data=cluster_count)
    data.tables[key] = t


def make_pipeline(target: Path):
    pipeline = import_file(target / "topo.data", atom_style="angle")

    # Load trajectory:
    mod = LoadTrajectoryModifier()
    pipeline.modifiers.append(mod)
    mod.source.load(str(target / "traj.nc"))
    return pipeline


def calc_hfield(target: Path, n_modes: int = 20, sel=slice(None, None, None)):

    def get_H(data, N):
        particle_cluster = np.array(data.particles["Cluster"])

        cell = data.cell
        boxsize = np.array([cell[i, i] for i in range(3)])
        boxorigin = np.array(cell[:, 3])

        particle_pos = np.array(data.particles["Position"])
        particle_sel = particle_cluster == 1
        pos = particle_pos[particle_sel]

        # shape = np.ceil(boxsize[:2] / delta).astype(int)
        shape = (N, N)

        bins = [
            np.linspace(boxorigin[i], boxorigin[i] + boxsize[i], num=1 + shape[i])
            for i in range(2)
        ]

        inds = [
            np.clip(np.digitize(pos[:, i], bins[i]), 1, len(bins[i]) - 1) - 1
            for i in range(2)
        ][::-1]
        IND = np.stack(inds, axis=1)

        pos_z = pos[:, 2]

        Zs = np.zeros(shape)
        # print(IND.shape)
        np.add.at(Zs, tuple(IND.T), pos_z)

        Cs = np.zeros(shape, dtype=int)
        np.add.at(Cs, tuple(IND.T), 1)

        assert (~(Cs < 10)).all()
        Z = Zs / Cs
        return Z, boxsize[0]

    pipeline = make_pipeline(target)

    pipeline.modifiers.append(ExpressionSelectionModifier(expression="ParticleType==2"))

    pipeline.modifiers.append(
        ClusterAnalysisModifier(cutoff=1.5, only_selected=True, sort_by_size=True)
    )
    pipeline.modifiers.append(ClearSelectionModifier())

    pipeline.modifiers.append(modify_clusters)

    frame_sel = np.arange(pipeline.num_frames)[sel]

    l = []
    l_step = []
    l_time = []
    for frame in tqdm(frame_sel):
        data = pipeline.compute(frame)
        l_step.append(data.attributes["Step"])
        l_time.append(data.attributes["Time"])
        l.append(get_H(data, n_modes))
    Z = np.stack([Z for Z, L in l], axis=0)
    L = np.stack([L for Z, L in l], axis=0)

    return xr.Dataset(
        {
            "z": (("step", "x", "y"), Z),
            "L": (("step",), L),
        },
        coords={"step": np.array(l_step), "time": ("step", np.array(l_time))},
    )


def thermo_2_df(target: Path, debug: bool = False):
    with target.open("rt") as f:
        _l = 0

        def read_thermo_run(o):
            nonlocal _l
            header = next(o).split()
            _l += 1
            # strip leading nulls
            r = []
            _l2 = []
            for l in o:
                _l += 1
                l = l.strip("\x00")
                if l.startswith("WARNING") or l.startswith("Fix halt"):
                    pass
                else:
                    if l.startswith("Loop time") or not l:
                        break
                    try:
                        le = [nan if e == "null" else float(e) for e in l.split()]
                        if not le:
                            break
                        r.append(le)
                        _l2.append(_l)
                    except:
                        break
            df = pd.DataFrame(data=r, columns=header)
            df["Step"] = df["Step"].astype(int)
            df = df.rename(columns={"Step": "step"})
            df = df.set_index("step")
            if debug:
                df["line"] = _l2
            return df

        def find_start(o):
            nonlocal _l
            for l in o:
                _l += 1
                if l.startswith("Per MPI rank memory"):
                    return True

            return False

        runs: list[pd.DataFrame] = []
        for _ in f:
            _l += 1
            if find_start(f):
                x = read_thermo_run(f)
                if debug:
                    x["run"] = len(runs)
                runs.append(x)

        if len(runs) == 0:
            columns = ["Time", "CPU", "step"]
            if debug:
                columns.extend(["line", "run"])
            df = pd.DataFrame(columns=columns)
        else:
            df = pd.concat(runs)
        df = df.reset_index()
        return df


def div0[T: np.ndarray](a: T, b: T, fill: float = np.nan) -> T:
    """a / b, divide by 0 -> `fill`
    div0( [-1, 0, 1], 0, fill=np.nan) -> [nan nan nan]
    div0( 1, 0, fill=np.inf ) -> inf
    """
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        c = np.true_divide(a, b)
    if np.isscalar(c):
        return c if np.isfinite(c) else fill  # type:ignore
    else:
        c[~np.isfinite(c)] = fill
        return c  # type:ignore


@dataclass
class mlh_res:
    popt: np.ndarray
    perr: np.ndarray
    Q: float

    @staticmethod
    def nan_from_p0(p0: np.ndarray):
        return mlh_res(
            popt=np.full_like(p0, np.nan), perr=np.full_like(p0, np.nan), Q=np.nan
        )


def mlh(
    f: Callable,
    p0: np.ndarray,
    bounds: tuple[np.ndarray, ...],
    x: np.ndarray,
    y: np.ndarray,
    e: np.ndarray,
):
    popt: np.ndarray
    pcov: np.ndarray
    e_valid = not (np.isnan(e).all())

    # drop points without errors
    if e_valid:
        sel = np.isfinite(e)
        x, y, e = x[sel], y[sel], e[sel]

    if len(x) == 0:
        return mlh_res.nan_from_p0(p0)
    try:
        popt, pcov = scipy.optimize.curve_fit(
            f,
            x,
            y,
            p0=p0,
            sigma=e if e_valid else None,
            absolute_sigma=e_valid,
            bounds=bounds,
        )[:2]
    except RuntimeError:
        return mlh_res.nan_from_p0(p0)
    pcov_dia: np.ndarray = np.diag(pcov)
    perr: np.ndarray = np.sqrt(pcov_dia)

    chisq = np.sum(div0((y - f(x, *popt)), e) ** 2)
    N = len(p0)
    M = len(x)
    k = M - N
    # print(N,M,k,chisq)

    Q: float = scipy.special.gammaincc(k / 2, chisq / 2) / scipy.special.gamma(k / 2)

    return mlh_res(popt, perr, Q)


class Model(Protocol):
    func: Callable
    p0: np.ndarray
    bounds: tuple

    def fit(self, x: np.ndarray, y: np.ndarray, e: np.ndarray) -> mlh_res: ...


def f_kat(x, k, a, t):
    return div0(1.0, k * x**a + t * x**2)


def _fit(self, x: np.ndarray, y: np.ndarray, e: np.ndarray):
    return mlh(f=self.func, p0=self.p0, bounds=self.bounds, x=x, y=y, e=e)


@dataclass
class Model_f_u_k(Model):
    func: Callable = lambda x, k: cast(float, div0(1, k * x**4, fill=np.nan))
    p0 = np.array(
        [
            10,
        ]
    )
    bounds = tuple(np.array([[0.0, np.inf]]).T)
    fit = _fit


def test_f_u_k():
    m = 10
    x = np.logspace(0.1, 1, base=10, num=m)
    p = (4,)
    m = Model_f_u_k()
    y = m.func(x, *p)
    e = m.func(x * (1 + np.random.default_rng().random(x.shape) * 1e-4), *p) - y
    r = m.fit(x=x, y=y, e=e)
    assert np.allclose(np.array([*p]), r.popt) and r.Q > 1e-2, r


@dataclass
class fit_Q_res:
    mlh_res: mlh_res
    x_last: float


def fit_Q(
    x: np.ndarray,
    y: np.ndarray,
    e: np.ndarray,
    m: Model,
    Q_threshold: float = 1e-4,
    x_max: float = np.inf,
):
    n = len(x)
    js = np.arange(n)[(x <= x_max)]
    js = js[js >= 3]  # at least 3 points must be present
    assert len(js) > 0, repr(js)

    fit_js = [m.fit(x=x[:j], y=y[:j], e=e[:j]) for j in js]

    Qs = np.array([fit.Q for fit in fit_js])
    Qgood = np.where(Qs > Q_threshold)[0]

    if len(Qgood) == 0:
        x_last = np.inf
        fit_res = mlh_res.nan_from_p0(m.p0)
    else:
        k = int(Qgood[-1])
        fit_res = fit_js[k]
        x_last = x[js[k]]
    return fit_Q_res(mlh_res=fit_res, x_last=x_last)


models = {
    "k": Model_f_u_k,
}

import matplotlib.pyplot as plt
from uncertainties import ufloat


def fs_fit_df(
    df: pd.DataFrame,
    do_plot: bool = False,
    plot_x_range=(10**-1, 10**0),
    fit_x_range=(2 * np.pi / 32, 2 * np.pi / 12),
):
    """_summary_

    Parameters
    ----------
    df : pd.DataFrame
        index -> q,v,e
    ...

    Returns
    -------

    """
    plot_fit_x_range = fit_x_range
    if do_plot:
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim(*plot_x_range)  # right = x_max if np.isfinite(x_max) else None)

        plt.axvline(x=fit_x_range[1], color="black", alpha=0.5)
        plt.axvline(x=fit_x_range[0], color="black", alpha=0.5)

    if "u" in df:
        x = np.array(df["q"], dtype=float)
        y, e = np.array([e.n for e in df["u"]], dtype=float), np.array(
            [e.s for e in df["u"]], dtype=float
        )
    else:
        x, y, e = df.reset_index()[["q", "v", "e"]].to_numpy().T

    if do_plot:
        plt.errorbar(
            x=x, y=y, yerr=e, ecolor="black", linestyle="none", marker="o", markersize=1
        )

    s = (fit_x_range[0] <= x) & (x <= fit_x_range[1])
    x, y, e = x[s], y[s], e[s]

    model_fit = {}

    m = Model_f_u_k()
    fit_res = m.fit(x=x, y=y, e=e)

    param_names = ["k"]
    model_fit["Q"] = fit_res.Q
    model_fit |= {
        k: ufloat(v, e) for k, v, e in zip(param_names, fit_res.popt, fit_res.perr)
    }

    if do_plot:
        x2 = np.logspace(np.log10(plot_fit_x_range[0]), np.log10(plot_fit_x_range[1]))
        plt.plot(
            x2,
            m.func(x2, *fit_res.popt),
            label="theory",
            color="grey",
            alpha=0.5,
            linestyle="--",
        )

    if do_plot:
        plt.legend()
        plt.xlabel(R"Wavenumber $q / \sigma^{-1}$")
        plt.ylabel(R"Scaled power $\langle h^2 \rangle L^2 / \sigma^4$")

    return model_fit


def main():
    # compute height field

    import sys

    target = Path(sys.argv[-1])

    # get average box size from thermo data
    df = thermo_2_df(target / "log.lammps")
    df["f_avg_lx"] = np.where(df["f_avg_lx"] == 0.0, np.nan, df["f_avg_lx"])

    step_discard = (
        10e3 / 60e3 * df.step.max()
    )  # proportionally remove initial trajectory for equilibration. for shorter than 60e3 it will not be accurate.

    plt.figure()
    plt.plot(df["step"] * 1e-2, df["f_avg_lx"])
    plt.xlabel(R"$\tau$")
    plt.ylabel(R"$L$")
    plt.axvline(
        x=step_discard * 1e-2,
        linestyle="--",
        alpha=0.5,
        label="equilibration threshold",
    )

    plt.legend()
    plt.savefig("L.pdf")
    plt.clf()

    lx_ufloat = series2ufloat(
        np.array(df[df.step > step_discard].iloc[:-1]["f_avg_lx"])
    )
    assert lx_ufloat is not None

    L_avg = lx_ufloat.n
    print(f"L_avg={lx_ufloat}")
    n_modes = 2 * (int(L_avg / 3) // 2)

    # compute height field
    Z_path = target / f"H_{n_modes}.nc"
    if not (Z_path).exists():
        Z = calc_hfield(target, n_modes=n_modes)
        Z.to_netcdf(Z_path)
    z_ds = (
        xr.load_dataset(Z_path)
        .sel(step=slice(step_discard, None))
        .swap_dims({"step": "time"})
    )
    # fourier transform it and average powers
    u = zL2u(z_ds, L_avg=L_avg)
    spectrum = (
        u.to_dataframe()
        .drop(columns=["nx", "ny"])
        .rename(columns={"nn": "n"})
        .reset_index()
        .sort_values("n")
    )
    spectrum.to_csv(target / "spectrum.csv")

    spectrum = spectrum.query("ok")

    # plot and fit
    model_fit = fs_fit_df(spectrum, do_plot=True)
    print("\n".join(f"{k}={v}" for k, v in model_fit.items()))

    plt.savefig("spectrum.pdf")
    # TODO: move code for multiple replica merging


if __name__ == "__main__":
    main()
