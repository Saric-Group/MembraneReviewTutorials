from typing import Hashable

import numpy as np
import xarray as xr
from equilibration import statisticalInefficiency_multiscale


def xr_norm(obj: xr.DataArray, dim: Hashable, ord=None):
    return xr.apply_ufunc(
        np.linalg.norm, obj, input_core_dims=[[dim]], kwargs={"ord": ord, "axis": -1}
    )


def xr_continuous_rfft2(
    v: xr.DataArray,
    input_dims=("x", "y"),
    output_dims=("nx", "ny"),
):
    """RFFT2 for xarray

    Parameters
    ----------
    v : xr.DataArray
        [x,y]->float
    input_dims : tuple, optional
        by default ("x", "y")
    output_dims : tuple, optional
        by default ("nx", "ny")

    Returns
    -------
    xr.DataArray
        [nx,ny]->float
        Fourier transform of v
    """

    lxly = (v.sizes[input_dims[0]], v.sizes[input_dims[1]])
    lx, ly = lxly
    assert lx % 2 == 0 and ly % 2 == 0  # TODO: fix
    nx = (np.fft.rfftfreq(lx, d=1) * lx).astype(int)
    ny = np.fft.fftshift(np.fft.fftfreq(ly, d=1) * ly).astype(int)

    r: xr.DataArray = xr.apply_ufunc(
        lambda x: np.fft.fftshift(
            np.fft.rfft2(x, axes=(-1, -2), norm="forward"), axes=-1
        ),
        v,
        input_core_dims=[[*input_dims]],
        output_core_dims=[[*output_dims]],
    )

    r = r.assign_coords(
        {output_dims[0]: (output_dims[0], nx), output_dims[1]: (output_dims[1], ny)}
    )

    return r


def xr_continuous_irfft2(
    v: xr.DataArray,
    s: tuple[int, int] | None = None,
    input_dims=("nx", "ny"),
    output_dims=("x", "y"),
):
    nx = v.sizes[input_dims[0]]
    ny = v.sizes[input_dims[1]]

    lxly = (2 * (nx - 1), ny) if s is None else s
    lx, ly = lxly

    v = v.reindex(
        {
            input_dims[0]: (np.fft.rfftfreq(lx, d=1.0) * lx).astype(int),
            input_dims[1]: np.fft.fftshift(np.fft.fftfreq(ly, d=1.0) * ly).astype(int),
        },
        fill_value=0,  # type:ignore
        copy=False,
    )

    r: xr.DataArray = xr.apply_ufunc(
        lambda x: np.fft.irfft2(
            np.fft.ifftshift(x, axes=-1), s=lxly[::-1], axes=(-1, -2), norm="forward"
        ),
        v,
        input_core_dims=[[*input_dims]],
        output_core_dims=[[*output_dims]],
    )
    return r


def xr_is_nyquist(
    Z: xr.DataArray, s: tuple[int, int], dims=("nx", "ny")
) -> xr.DataArray:
    lx, ly = s
    ca, cb = (lx % 2 == 0), (ly % 2 == 0)

    a = Z[dims[0]] == lx / 2
    b = Z[dims[1]] == -(ly / 2)

    if ca and cb:
        return a | b
    elif ca:
        return a
    elif cb:
        return b
    else:
        return xr.DataArray(False)


def calc_spectrum_grid_shift(Z: xr.DataArray, r: xr.DataArray) -> xr.DataArray:
    nx, ny = Z["nx"], Z["ny"]
    n = xr.concat([nx, ny], "spatial")
    k = 2 * np.pi * n
    a = xr.DataArray(
        [1 / r.sizes["x"], 1 / r.sizes["y"]],
        dims=("spatial",),
        coords={"spatial": ["x", "y"]},
    )
    return np.exp(0.5j * xr.dot(k, a, dims="spatial"))


def xr_spectrum_from_grid(z: xr.DataArray) -> xr.DataArray:
    """
    Computes Fourier transform of a grid-sampled function in [0,1]^2
    Input:  xr.DataArray(values,coords={'x','y'}), coords being regularly spaced on each dimension
    Output: xr.DataArray(complex values, coords = {'nx','ny'}) ,'n*' being the wavenumbers
    """

    Z = xr_continuous_rfft2(z)

    Z /= calc_spectrum_grid_shift(Z, z)

    return Z


def test_xr_spectrum_from_grid():
    def make_z(
        L: float = 1.0, N: int = 2**5, nx: int = 1, ny: int = 1, amp: float = 1.0
    ) -> xr.DataArray:
        ls = np.linspace(0, L, num=N + 1)
        ps = (ls[1:] + ls[:-1]) / 2

        # centers of 2D N by N regular lattice in [0,1]^2
        r = xr.DataArray(
            np.stack(np.meshgrid(ps, ps, indexing="ij")), dims=("spatial", "x", "y")
        ).assign_coords(x=ps, y=ps, spatial=["x", "y"])

        # nx, ny = 2, 1
        n = xr.DataArray([nx, ny], coords={"spatial": ["x", "y"]}, dims=("spatial",))
        k = 2 * np.pi / (L / n)

        z = amp * xr.apply_ufunc(np.cos, (k * r).sum("spatial"))

        return z

    L = 1.0
    N = 2**5
    nx = 1
    ny = 0
    amp = 1.0
    z = make_z(L=L, N=N, nx=nx, ny=ny, amp=amp)

    Z = xr_spectrum_from_grid(z)
    Z = Z.where(~xr.apply_ufunc(np.isclose, Z, 0), drop=True)
    assert Z.shape == (1, 1)
    assert Z.nx[0] == nx and Z.ny[0] == ny
    assert np.allclose(np.abs(Z), 1 / 2)
    assert np.allclose(np.angle(Z), 0)


def Z2g(Z: xr.DataArray):
    """
    Extracts equilibrated wave powers from complex spectrum, by removing initial non-stable sequence and rescaling std. dev. according to correlation time, taking in account phase and amplitude, for each wavenumber.
    TODO: autocorrelation of complex time series?
    """

    output = {
        "v": float,
        "e": float,
        "g": float,
    }
    output_names = tuple(output.keys())
    output_dtypes = list(output.values())

    def _np_measure(Z):

        na_res = np.nan, np.nan, 0, False, 0, np.nan
        if np.isnan(Z).any():
            return na_res
        Zeng = np.abs(Z) ** 2
        Zang = np.angle(Z)

        def method(x):
            ineff = statisticalInefficiency_multiscale(x)
            if ineff is None:
                return None
            elif np.isinf(ineff):
                return None
            else:
                return ineff

        ineffs = [method(x) for x in (Zeng, Zang)]
        n = len(Z)
        ineff = None if any(i is None for i in ineffs) else max(ineffs)
        if ineff is None:
            return na_res
        zs = Zeng
        return zs.mean(), (zs.std() / (len(zs) - 1) ** 0.5) * (ineff) ** 0.5, ineff

    g = xr.Dataset(
        dict(
            zip(
                output_names,
                xr.apply_ufunc(
                    _np_measure,
                    Z,
                    input_core_dims=[["time"]],
                    output_core_dims=tuple(() for _ in output_names),
                    output_dtypes=output_dtypes,
                    vectorize=True,
                ),
            )
        )
    )
    return g


def z2Z(z: xr.DataArray):
    """Computes the amplitude squared vs wavenumber from height field timeseries in the form of a rectilinear grid obtained by averaging height in each cell.
    Cell lengths are not used; scale

    Parameters
    ----------
    z : xr.DataArray
        [time,x,y] -> float
    Returns
    -------
    u : xr.Dataset
        Dimensions
            n : float
                wavenumber
        Variables
            v : [n] : float
                amplitude squared (A^2)
            e : [n] : float
                std dev of the mean of A^2
    """

    Z = xr_spectrum_from_grid(z)
    Z = Z.assign_coords(
        nn=xr_norm(xr.concat([Z["nx"], Z["ny"]], "spatial").rename("n"), "spatial")
    )
    assert not np.isnan(Z).any()
    Z = (
        Z.where((Z.nn > 0) & (~xr_is_nyquist(Z, (z.sizes["x"], z.sizes["y"]))))
        .stack(nxy=["nx", "ny"])
        .dropna("nxy")
    )
    assert not np.isnan(Z).any()

    # correct for averaging
    Z /= xr.apply_ufunc(
        np.sinc,
        Z["nx"] / z.sizes["x"],
    ) * xr.apply_ufunc(
        np.sinc,
        Z["ny"] / z.sizes["y"],
    )

    return Z


def zL2u(ds: xr.Dataset, L_avg: float):
    """Computes ( A(amplitude)*avg(L(box length)) )^2 vs wavenumber from height field timeseries in the form of a rectilinear grid obtained by averaging height in each cell.
    Groups measurements by wavenumber, and assumes box has approximately constant length.

    Parameters
    ----------
    z : xr.DataArray
        [time,x,y] -> float
        timeseries of 2d height field
    L_avg : float
        average box size

    Returns
    -------
    u : xr.Dataset
        Dimensions
            n : float
                wavenumber
        Variables
            v : [n] : float
                average over time of amplitude squared
            e : [n] : float
                standard deviation of the mean of the amplitude squared
    """

    z = ds["z"]
    # fourier transform z_{i,j} to Z_{i,j}
    Z = z2Z(z)

    # only trajectories with constant spacing between samples
    assert ((d := np.diff(Z.time.data)) == d[0]).all()
    delta_time = d[0]

    ZL = Z * L_avg

    u = Z2g(ZL)
    u["tau"] = (u["g"] - 1) * delta_time / 2
    u["ok"] = (Z.sizes["time"] / u["g"]) > 20
    u["q"] = 2 * np.pi * u["nn"] / L_avg
    return u
