use core::f32;

use crate::ray::ray::{Hit, Ray, RayCast};

pub struct Sphere {
    pub cnt: [f32; 3],
    pub r: f32,
}

impl RayCast for Sphere {
    fn raycast(&self, ray: &Ray) -> Option<Hit> {
        if let Some(t) = del_geo_nalgebra::sphere::intersection_ray(
            &nalgebra::Vector3::<f32>::from(self.cnt),
            self.r,
            &nalgebra::Vector3::<f32>::from(ray.o),
            &nalgebra::Vector3::<f32>::from(ray.d),
        ) {
            return Some(Hit {
                ray: ray.clone(),
                t,
            });
        }

        None
    }
}

/// return hit, index of sphere
pub fn ray_cast_shperes(
    (ix, iy): (usize, usize),
    (dx, dy): (f32, f32),
    img_shape: (usize, usize),
    fov: f32,
    transform_camlcl2world: &[f32; 16],
    spheres: &[Sphere],
) -> Option<(Hit, usize)> {
    use crate::cam;
    let (ray_org, ray_dir) =
        cam::gen_ray((ix, iy), (dx, dy), img_shape, fov, transform_camlcl2world);
    let ray = Ray {
        o: ray_org,
        d: ray_dir,
    };

    let (hit, isphere) = {
        let mut hit = Hit {
            ray: ray.clone(),
            t: f32::INFINITY,
        };
        let mut isphere = 0;
        for (i, sphere) in spheres.iter().enumerate() {
            if let Some(h) = sphere.raycast(&ray) {
                if h.t < hit.t {
                    hit.t = h.t;
                    isphere = i;
                }
            }
        }

        (hit, isphere)
    };

    if hit.t == f32::INFINITY {
        return None;
    }

    Some((hit, isphere))
}
