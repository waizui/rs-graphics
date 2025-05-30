use num_traits::Float;

pub fn lerp<T>(a: &[T; 3], b: &[T; 3], t: T) -> [T; 3]
where
    T: Float,
{
    let v0 = (T::one() - t) * a[0] + (t) * b[0];
    let v1 = (T::one() - t) * a[1] + (t) * b[1];
    let v2 = (T::one() - t) * a[2] + (t) * b[2];
    [v0, v1, v2]

}
