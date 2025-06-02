#[derive(Clone)]
pub struct Ray {
    pub o: [f32; 3],
    pub d: [f32; 3],
}

pub struct Hit {
    pub ray: Ray,
    pub t: f32,
}

pub trait RayCast {
    fn raycast(&self, ray: &Ray) -> Option<Hit>;
}
