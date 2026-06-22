use super::*;

// ---- Web Mercator (spherical) ----

/// Backend-typed Web Mercator constants.
///
/// The only spherical Mercator projection constant is the sphere radius. This
/// small wrapper keeps the API consistent with the other projection families:
/// lift or construct the backend radius once, then reuse it across a batch.
#[derive(Clone, Debug)]
pub struct WebMercatorConstants<T> {
    /// Spherical radius in metres.
    pub radius: T,
}

impl<T: Real> WebMercatorConstants<T> {
    /// Build constants from a backend-typed radius.
    #[must_use]
    pub fn from_radius(radius: T) -> Self {
        Self { radius }
    }
}

/// Web Mercator (EPSG:3857) forward: geographic lat/lon **degrees** to projected
/// `(x, y)` metres on a sphere of `radius` metres.
///
/// `y = R·ψ(φ)` with `ψ` the isometric latitude `log((1 + sinφ)/cosφ)` — the
/// exact spherical-Mercator form (`riemann` `Hausdorff/Rhumb.isoLat`).
#[must_use]
pub fn web_mercator_forward<T: Real>(lat_deg: T, lon_deg: T, radius: T) -> (T, T) {
    let c = WebMercatorConstants::from_radius(radius);
    web_mercator_forward_typed(lat_deg, lon_deg, &c)
}

/// Web Mercator forward with already-typed constants.
#[must_use]
pub fn web_mercator_forward_typed<T: Real>(
    lat_deg: T,
    lon_deg: T,
    c: &WebMercatorConstants<T>,
) -> (T, T) {
    let one = T::one();
    let (s, cos_phi) = lat_deg.to_radians().sin_cos();
    let x = c.radius.clone() * lon_deg.to_radians();
    let y = c.radius.clone() * ((one + s) / cos_phi).ln();
    (x, y)
}

/// Web Mercator inverse: projected `(x, y)` metres back to geographic lat/lon
/// **degrees**, on a sphere of `radius` metres.
///
/// `φ = arcsin(tanh(y/R))` — the provably exact inverse (`riemann`
/// `Hausdorff/Mercator.mercator_inverse`), since `tanh(ψ(φ)) = sinφ`.
#[must_use]
pub fn web_mercator_inverse<T: Real>(x: T, y: T, radius: T) -> (T, T) {
    let c = WebMercatorConstants::from_radius(radius);
    web_mercator_inverse_typed(x, y, &c)
}

/// Web Mercator inverse with already-typed constants.
#[must_use]
pub fn web_mercator_inverse_typed<T: Real>(x: T, y: T, c: &WebMercatorConstants<T>) -> (T, T) {
    let lat_deg = (y / c.radius.clone()).tanh().asin().to_degrees();
    let lon_deg = (x / c.radius.clone()).to_degrees();
    (lat_deg, lon_deg)
}

/// Project a [`LatLon`] provider through [`web_mercator_forward`] — the entry
/// point a tensor/columnar lat/lon source uses to project en masse.
#[must_use]
pub fn web_mercator_project_at<T, P>(p: &P, radius: T) -> (T, T)
where
    T: Real,
    P: LatLon<T>,
{
    web_mercator_forward(p.lat_deg(), p.lon_deg(), radius)
}

/// Project a [`LatLon`] provider through [`web_mercator_forward_typed`].
#[must_use]
pub fn web_mercator_project_at_typed<T, P>(p: &P, c: &WebMercatorConstants<T>) -> (T, T)
where
    T: Real,
    P: LatLon<T>,
{
    web_mercator_forward_typed(p.lat_deg(), p.lon_deg(), c)
}
