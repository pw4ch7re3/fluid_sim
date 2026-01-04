macro_rules! ix {
    ($x:expr, $y:expr, $z: expr, $n: expr) => {
        $x as usize +
        $y as usize * $n as usize +
        $z as usize * $n as usize * $n as usize
    };
}

#[allow(unused)]
struct FluidCube {
    size: i32,
    dt: f32,
    diff: f32,
    visc: f32,

    s: Vec<f32>,
    density: Vec<f32>,

    vx: Vec<f32>,
    vy: Vec<f32>,
    vz: Vec<f32>,

    vx0: Vec<f32>,
    vy0: Vec<f32>,
    vz0: Vec<f32>,
}

#[allow(unused)]
impl FluidCube {
    fn new(size: i32, diffusion: f32,
           viscosity: f32, dt: f32) -> FluidCube
    {
        let n = size as usize;
        let dim = n * n * n;
        FluidCube {
            size,
            dt,
            diff: diffusion,
            visc: viscosity,
            
            s: vec![0.0; dim],
            density: vec![0.0; dim],
            
            vx: vec![0.0; dim],
            vy: vec![0.0; dim],
            vz: vec![0.0; dim],
            
            vx0: vec![0.0; dim],
            vy0: vec![0.0; dim],
            vz0: vec![0.0; dim],
        }
    }

    fn add_density(&mut self, x: i32, y: i32, z: i32, amount: f32) {
        self.density[ix!(x, y, z, self.size)] += amount;
    }
    
    fn add_velocity(&mut self, x: i32, y: i32, z: i32,
                    amount_x: f32, amount_y: f32, amount_z: f32)
    {
        let index = ix!(x, y, z, self.size);
        self.vx[index] += amount_x;
        self.vy[index] += amount_y;
        self.vz[index] += amount_z;
    }

    fn set_bnd(b: i32, x: &mut Vec<f32>, n: i32) {
        for j in 1..(n - 1) {
            for i in 1..(n - 1) {
                x[ix!(i, j, 0, n)] =
                    if b == 3 { -x[ix!(i, j, 1, n)] }
                    else { x[ix!(i, j, 1, n)] };
                x[ix!(i, j, n - 1, n)] =
                    if b == 3 { -x[ix!(i, j, n - 2, n)] }
                    else { x[ix!(i, j, n - 2, n)] };
            }
        }
        for k in 1..(n - 1) {
            for i in 1..(n - 1) {
                x[ix!(i, 0, k, n)] =
                    if b == 2 { -x[ix!(i, 1, k, n)] }
                    else { x[ix!(i, 1, k, n)] };
                x[ix!(i, n - 1, k, n)] =
                    if b == 2 { -x[ix!(i, n - 2, k, n)] }
                    else { x[ix!(i, n - 2, k, n)] };
            }
        }
        for k in 1..(n - 1) {
            for j in 1..(n - 1) {
                x[ix!(0, j, k, n)] =
                    if b == 1 { -x[ix!(1, j, k, n)] }
                    else { x[ix!(1, j, k, n)] };
                x[ix!(n - 1, j, k, n)] =
                    if b == 1 { -x[ix!(n - 2, j, k, n)] }
                    else { x[ix!(n - 2, j, k, n)] };
            }
        }
    
        const ONETHIRD: f32 = 0.33;
        x[ix!(0, 0, 0, n)] =
            ONETHIRD * (x[ix!(1, 0, 0, n)] +
                        x[ix!(0, 1, 0, n)] +
                        x[ix!(0, 0, 1, n)]);
        x[ix!(n - 1, 0, 0, n)] =
            ONETHIRD * (x[ix!(n - 2, 0, 0, n)] +
                        x[ix!(n - 1, 1, 0, n)] +
                        x[ix!(n - 1, 0, 1, n)]);
        x[ix!(0, n - 1, 0, n)] =
            ONETHIRD * (x[ix!(1, n - 1, 0, n)] +
                        x[ix!(0, n - 2, 0, n)] +
                        x[ix!(0, n - 1, 1, n)]);
        x[ix!(0, 0, n - 1, n)] =
            ONETHIRD * (x[ix!(1, 0, n - 1, n)] +
                        x[ix!(0, 1, n - 1, n)] +
                        x[ix!(0, 0, n - 2, n)]);
        x[ix!(n - 1, n - 1, 0, n)] =
            ONETHIRD * (x[ix!(n - 2, n - 1, 0, n)] +
                        x[ix!(n - 1, n - 2, 0, n)] +
                        x[ix!(n - 1, n - 1, 1, n)]);
        x[ix!(n - 1, 0, n - 1, n)] =
            ONETHIRD * (x[ix!(n - 2, 0, n - 1, n)] +
                        x[ix!(n - 1, 1, n - 1, n)] +
                        x[ix!(n - 1, 0, n - 2, n)]);
        x[ix!(0, n - 1, n - 1, n)] =
            ONETHIRD * (x[ix!(1, n - 1, n - 1, n)] +
                        x[ix!(0, n - 2, n - 1, n)] +
                        x[ix!(0, n - 1, n - 2, n)]);    
        x[ix!(n - 1, n - 1, n - 1, n)] =
            ONETHIRD * (x[ix!(n - 2, n - 1, n - 1, n)] +
                        x[ix!(n - 1, n - 2, n - 1, n)] +
                        x[ix!(n - 1, n - 1, n - 2, n)]);
    }
    
    fn lin_solve(b: i32, x: &mut Vec<f32>, x0: &Vec<f32>,
                 a: f32, c: f32, iter: i32, n: i32)
    {
        let recip_c = 1.0 / c;
        for _ in 0..iter {
            for k in 1..(n - 1) {
                for j in 1..(n - 1) {
                    for i in 1..(n - 1) {
                        x[ix!(i, j, k, n)] = recip_c * (
                            x0[ix!(i, j, k, n)] +
                            a * (x[ix!(i + 1, j, k, n)] +
                                 x[ix!(i - 1, j, k, n)] +
                                 x[ix!(i, j + 1, k, n)] +
                                 x[ix!(i, j - 1, k, n)] +
                                 x[ix!(i, j, k + 1, n)] +
                                 x[ix!(i, j, k - 1, n)])
                        );
                    }
                }
            }
            Self::set_bnd(b, x, n);
        }
    }
    
    fn diffuse(b: i32, x: &mut Vec<f32>, x0: &Vec<f32>,
               diff: f32, dt: f32, iter: i32, n: i32)
    {
        let a: f32 = dt * diff * (n - 2) as f32 * (n - 2) as f32;
        Self::lin_solve(b, x, x0, a, (1.0 + 6.0 * a), iter, n);
    }
    
    fn project(veloc_x: &mut Vec<f32>,
               veloc_y: &mut Vec<f32>,
               veloc_z: &mut Vec<f32>,
               p: &mut Vec<f32>, div: &mut Vec<f32>, iter: i32, n: i32)
    {
        for k in 1..(n - 1) {
            for j in 1..(n - 1) {
                for i in 1..(n - 1) {
                    div[ix!(i, j, k, n)] = -0.5f32 * (
                        veloc_x[ix!(i + 1, j, k, n)] -
                        veloc_x[ix!(i - 1, j, k, n)] +
                        veloc_y[ix!(i, j + 1, k, n)] -
                        veloc_y[ix!(i, j - 1, k, n)] +
                        veloc_z[ix!(i, j, k + 1, n)] -
                        veloc_z[ix!(i, j, k - 1, n)]
                    ) / n as f32;
                    p[ix!(i, j, k, n)] = 0.0;
                }
            }
        }
        Self::set_bnd(0, div, n);
        Self::set_bnd(0, p, n);
        Self::lin_solve(0, p, div, 1.0, 6.0, iter, n);
    
        for k in 1..(n - 1) {
            for j in 1..(n - 1) {
                for i in 1..(n - 1) {
                    veloc_x[ix!(i, j, k, n)] -= 0.5f32 *
                        (p[ix!(i + 1, j, k, n)] -
                         p[ix!(i - 1, j, k, n)]) * n as f32;
                    veloc_y[ix!(i, j, k, n)] -= 0.5f32 *
                        (p[ix!(i, j + 1, k, n)] -
                         p[ix!(i, j - 1, k, n)]) * n as f32;
                    veloc_z[ix!(i, j, k, n)] -= 0.5f32 *
                        (p[ix!(i, j, k + 1, n)] -
                         p[ix!(i, j, k - 1, n)]) * n as f32;
                }
            }
        }
        Self::set_bnd(1, veloc_x, n);
        Self::set_bnd(2, veloc_y, n);
        Self::set_bnd(3, veloc_z, n);
    }
    
    fn advect(b: i32, d: &mut Vec<f32>, d0: &Vec<f32>,
              veloc_x: &Vec<f32>, veloc_y: &Vec<f32>, veloc_z: &Vec<f32>,
              dt: f32, n: i32)
    {
        let dt_x = dt * (n as f32 - 2.0);
        let dt_y = dt * (n as f32 - 2.0);
        let dt_z = dt * (n as f32 - 2.0);
    
        for k in 1..(n - 1) {
            for j in 1..(n - 1) {
                for i in 1..(n - 1) {
                    let dx = veloc_x[ix!(i, j, k, n)] * dt_x;
                    let dy = veloc_y[ix!(i, j, k, n)] * dt_y;
                    let dz = veloc_z[ix!(i, j, k, n)] * dt_z;
                    let x = i as f32 - dx;
                    let y = j as f32 - dy;
                    let z = k as f32 - dz;
    
                    let x = x.clamp(0.5f32,
                        n as f32 + 0.5f32);
                    let y = y.clamp(0.5f32,
                        n as f32 + 0.5f32);
                    let z = z.clamp(0.5f32,
                        n as f32 + 0.5f32);
    
                    let i0 = x.floor(); let i1 = i0 + 1.0f32;
                    let j0 = y.floor(); let j1 = j0 + 1.0f32;
                    let k0 = z.floor(); let k1 = k0 + 1.0f32;

                    let s1 = x - i0; let s0 = 1.0f32 - s1;
                    let t1 = y - j0; let t0 = 1.0f32 - t1;
                    let u1 = z - k0; let u0 = 1.0f32 - u1;
    
                    d[ix!(i, j, k, n)] =
                        s0 * (t0 * (u0 * d0[ix!(i0, j0, k0, n)] +
                                    u1 * d0[ix!(i0, j0, k1, n)]) +
                              t1 * (u0 * d0[ix!(i0, j1, k0, n)] +
                                    u1 * d0[ix!(i0, j1, k1, n)])) +
                        s1 * (t0 * (u0 * d0[ix!(i1, j0, k0, n)] +
                                    u1 * d0[ix!(i1, j0, k1, n)]) +
                              t1 * (u0 * d0[ix!(i1, j1, k0, n)] +
                                    u1 * d0[ix!(i1, j1, k1, n)]));
                }
            }
        }
        Self::set_bnd(b, d, n);
    }

    fn step(&mut self) {
        let n = self.size;
        let visc = self.visc;
        let diff = self.diff;
        let dt = self.dt;

        Self::diffuse(
            1, &mut self.vx, &self.vx0,
            visc, dt, 4, n
        );
        Self::diffuse(
            2, &mut self.vy, &self.vy0,
            visc, dt, 4, n
        );
        Self::diffuse(
            3, &mut self.vz, &self.vz0,
            visc, dt, 4, n
        );

        Self::project(
            &mut self.vx0,
            &mut self.vy0,
            &mut self.vz0,
            &mut self.vx, &mut self.vy, 4, n
        );

        Self::advect(
            1, &mut self.vx, &self.vx0,
            &self.vx0,
            &self.vy0,
            &self.vz0,
            dt, n
        );
        Self::advect(
            2, &mut self.vy, &self.vy0,
            &self.vx0,
            &self.vy0,
            &self.vz0,
            dt, n
        );
        Self::advect(
            3, &mut self.vz, &self.vz0,
            &self.vx0,
            &self.vy0,
            &self.vz0,
            dt, n
        );

        Self::project(
            &mut self.vx,
            &mut self.vy,
            &mut self.vz,
            &mut self.vx0, &mut self.vy0, 4, n
        );

        Self::diffuse(
            0, &mut self.s, &self.density,
            diff, dt, 4, n
        );
        Self::advect(
            0, &mut self.density, &self.s,
            &self.vx,
            &self.vy,
            &self.vz,
            dt, n
        );
    }
}

fn main() {
    let size = 16;
    let center = size / 2;

    let mut cube = FluidCube::new(
        size, 0.1, 0.0, 0.1
    );

    for i in 0..30 {
        cube.add_density(center, center, center, 100.0);
        cube.add_velocity(
            center, center, center,
            0.0, 2.0, 0.0
        );

        cube.step();

        let m: f32 = cube.density.iter().sum();
        let rho =
            cube.density[ix!(center, center, center, size)];
        println!("Step {i:02} | {m:>12.6} | {rho:>10.6}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let size = 10;
        let cube = FluidCube::new(
            size, 0.1, 0.1, 0.1
        );
        assert_eq!(cube.density.len(), (size * size * size) as usize);
        assert_eq!(cube.density[0], 0.0);
    }

    #[test]
    fn test_add_density() {
        let mut cube = FluidCube::new(
            10, 0.1, 0.1, 0.1
        );
        cube.add_density(5, 5, 5, 50.0);
        assert_eq!(cube.density[ix!(5, 5, 5, 10)], 50.0);
    }

    #[test]
    fn test_step() {
        let mut cube = FluidCube::new(
            8, 0.1, 0.1, 0.1
        );
        cube.add_density(4, 4, 4, 10.0);
        cube.step();
    }
}
