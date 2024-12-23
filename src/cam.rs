use del_geo_core::vec3;
use nalgebra::Matrix4;

type Real = f32;

#[derive(Debug)]
struct Camera {
    fov: Real,
    n: Real,
    f: Real,
    pos: [Real; 3],
    up: [Real; 3],
    forward: [Real; 3],
    matrix_v: Matrix4<Real>,
    matrix_p: Matrix4<Real>,
    matrix_vp: Matrix4<Real>,
}

impl Camera {
    pub fn new() -> Camera {
        todo!()
    }

    pub fn lookat(&mut self, forward: &[Real; 3]) {
        let mut up = [0., 1., 0.];
        let right = vec3::cross(&up, forward);
        up = vec3::cross(forward, &right);

        // formatter will not work if use matrix-style initializing
        let mut vmat = Matrix4::new(
            right[0], up[0], forward[0], 0., right[1], up[1], forward[1], 0., right[2], up[2],
            forward[2], 0., 0., 0., 0., 1.,
        )
        .transpose();

        vmat[(0, 3)] = -self.pos[0];
        vmat[(1, 3)] = -self.pos[1];
        vmat[(2, 3)] = -self.pos[2];

        self.up = up;
        self.forward = *forward;

        self.matrix_v = vmat;
        self.matrix_vp = self.matrix_p * vmat
    }
}

pub fn project_matrix(n: Real, f: Real, fov: Real) -> Matrix4<Real> {
    #[rustfmt::skip]
    let p = Matrix4::new(
        1., 0., 0., 0., 
        0., 1., 0., 0., 
        0., 0., f/(f-n), -f*n/(f-n), 
        0., 0., 1., 0.,
    );

    let invfov = 1. / (fov.to_radians() / 2.).tan();
    #[rustfmt::skip]
    let s = Matrix4::new(
        invfov, 0., 0., 0., 
        0., invfov, 0., 0., 
        0., 0., 1.,0. , 
        0., 0., 0., 1.,
    );

    s * p
}
