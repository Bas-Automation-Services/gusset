//! A cone or cylinder segment made on a press brake: planar facets between
//! the two arcs, joined by bends. Laid out flat as a chain of bends from the
//! first facet, which the fold places where it lies on the cone.

use opengeometry::brep::BrepEnvelope;

use crate::flat::{Cell, Curve, Cut, FlatBend, FlatPattern, Fold, Layout, Zone, P2};
use crate::fold::{Mass, Model, Rigid};
use crate::spec::{Cone, Facets, PartSpec, Radial};
use crate::SheetError;

type P3 = [f64; 3];

const EPS: f64 = 1e-9;

pub(crate) fn layout(spec: &PartSpec) -> Option<Result<Layout, SheetError>> {
    faceted(spec).map(|r| r.map(|(layout, _)| layout))
}

pub(crate) fn fold(spec: &PartSpec) -> Option<Result<BrepEnvelope, SheetError>> {
    faceted(spec).map(|r| {
        let (layout, place) = r?;
        crate::fold::solid(&Model::placed(spec, &layout, place))
    })
}

pub(crate) fn mass(spec: &PartSpec) -> Option<Result<Mass, SheetError>> {
    faceted(spec).map(|r| {
        let (layout, place) = r?;
        Ok(Model::placed(spec, &layout, place).exact_mass())
    })
}

fn faceted(spec: &PartSpec) -> Option<Result<(Layout, Rigid), SheetError>> {
    let cone = spec.base.cone.as_ref()?;
    let facets = cone.facets?;
    Some(build(spec, cone, facets))
}

fn build(spec: &PartSpec, cone: &Cone, facets: Facets) -> Result<(Layout, Rigid), SheetError> {
    crate::cone::check(spec, cone)?;
    let b = &spec.base;
    if !b.holes.is_empty() || !b.cutouts.is_empty() || !b.punches.is_empty() {
        return Err(SheetError::Unsupported(
            "holes, cut-outs and punches on a faceted cone".into(),
        ));
    }
    let (t, r, k) = (spec.thickness, spec.bend_radius, spec.k_factor);
    let [r0, r1] = cone.radii;
    let (h, sweep) = (cone.height, cone.sweep.to_radians());
    let n = facets
        .count(sweep, r0.max(r1))
        .ok_or_else(|| SheetError::Invalid("cone: the facets' limit must be positive".into()))?;
    let corner = |i: usize, radius: f64, z: f64| {
        let a = sweep * i as f64 / n as f64;
        [radius * a.cos(), radius * a.sin(), z]
    };
    let (p, q) = (|i| corner(i, r0, 0.0), |i| corner(i, r1, h));
    // Each facet's normal towards the axis: the bends' inside.
    let inward = |j: usize| unit(cross(sub(q(j), p(j)), sub(p(j + 1), p(j))));
    let angle = if n > 1 {
        let (a, c) = (inward(0), inward(1));
        norm(cross(a, c)).atan2(dot(a, c))
    } else {
        0.0
    };

    // The flat frame: +x the way the sheet unrolls, +y up the first ruling,
    // +z away from the axis; the arcs' face at z = z_s.
    let z_s = match cone.material {
        Radial::Outside => 0.0,
        Radial::Inside => t,
    };
    let y = unit(sub(q(0), p(0)));
    let z = scale(inward(0), -1.0);
    let x = cross(y, z);
    let place = Rigid::frame(sub(p(0), scale(z, z_s)), x, y, z);
    // A bend starts where its face on the arcs' side leaves the facet.
    let setback = (r + z_s) * (angle / 2.0).tan();
    let allowance = angle * (r + k * t);
    // Per facet: its plane's points to the flat, origin and axes.
    let mut develop = vec![(p(0), x, y, [0.0, 0.0], [1.0, 0.0], [0.0, 1.0])];
    let mut folds: Vec<Fold> = Vec::new();
    for j in 0..n - 1 {
        let (o, ex, ey, fo, fx, fy) = develop[j];
        let to_flat = |a: P3| {
            let d = sub(a, o);
            let (u, v) = (dot(d, ex), dot(d, ey));
            [fo[0] + u * fx[0] + v * fy[0], fo[1] + u * fx[1] + v * fy[1]]
        };
        let edge = to_flat(p(j + 1));
        let u = unit2(sub2(to_flat(q(j + 1)), edge));
        let normal = [u[1], -u[0]];
        folds.push(Fold {
            flange: format!("bend{}", j + 1),
            parent: j.checked_sub(1),
            edge_point: edge,
            u,
            n: normal,
            s0: -setback,
            bend_allowance: allowance,
            angle,
            radius: r,
            k_factor: k,
            up: false,
        });
        // The next facet: along the ruling and square to it, past the bend.
        let e1 = unit(sub(q(j + 1), p(j + 1)));
        let w = sub(p(j + 2), p(j + 1));
        let e2 = unit(sub(w, scale(e1, dot(w, e1))));
        let shift = allowance - 2.0 * setback;
        develop.push((
            p(j + 1),
            e1,
            e2,
            [edge[0] + shift * normal[0], edge[1] + shift * normal[1]],
            u,
            normal,
        ));
    }

    // Each bend's ends: where its face on the arcs' side meets z = 0 and
    // z = h, along the bend at flat s.
    let mut layout = Layout {
        flat: FlatPattern {
            outline: Vec::new(),
            outline_arcs: Vec::new(),
            holes: Vec::new(),
            cutouts: Vec::new(),
            bends: Vec::new(),
            punches: Vec::new(),
            forms: Vec::new(),
        },
        cells: Vec::new(),
        folds,
        holes: Vec::new(),
        embosses: Vec::new(),
    };
    let model = Model::placed(spec, &layout, place);
    let cuts: Vec<[Cut; 2]> = layout
        .folds
        .iter()
        .enumerate()
        .map(|(j, f)| {
            let rn = r + k * t;
            [0.0, h].map(|level| {
                let height = |s: f64, along: f64| {
                    let at = add2(f.edge_point, add2(scale2(f.n, s), scale2(f.u, along)));
                    model.map(Zone::Bend(j), at, z_s)[2]
                };
                let rows = [0.0, 0.5, 1.0].map(|part| {
                    let s = f.s0 + part * allowance;
                    let (z0, z1) = (height(s, 0.0), height(s, 1.0));
                    let (sin, cos) = ((s - f.s0) / rn).sin_cos();
                    [1.0, cos, sin, (level - z0) / (z1 - z0)]
                });
                Cut {
                    origin: f.edge_point,
                    u: f.u,
                    n: f.n,
                    s0: f.s0,
                    radius: rn,
                    c: crate::flat::solve3(rows),
                }
            })
        })
        .collect();

    // Facet j's corners: left and right, bottom and top.
    let corners = |j: usize| -> [P2; 4] {
        let (o, ex, ey, fo, fx, fy) = develop[j];
        let to_flat = |a: P3| {
            let d = sub(a, o);
            let (u, v) = (dot(d, ex), dot(d, ey));
            [fo[0] + u * fx[0] + v * fy[0], fo[1] + u * fx[1] + v * fy[1]]
        };
        let [lb, lt] = if j == 0 {
            [to_flat(p(0)), to_flat(q(0))]
        } else {
            cuts[j - 1].map(|c| c.point(c.s0 + allowance))
        };
        let [rb, rt] = if j == n - 1 {
            [to_flat(p(n)), to_flat(q(n))]
        } else {
            cuts[j].map(|c| c.point(c.s0))
        };
        [lb, rb, rt, lt]
    };
    // Each facet's chords must outlast its bends' setbacks.
    for (j, &(_, ex, ey, _, fx, fy)) in develop.iter().enumerate() {
        let [lb, rb, rt, lt] = corners(j);
        let along = sub(p(j + 1), p(j));
        let d = add2(scale2(fx, dot(along, ex)), scale2(fy, dot(along, ey)));
        if dot2(sub2(rb, lb), d) <= EPS || dot2(sub2(rt, lt), d) <= EPS {
            return Err(SheetError::Invalid(
                "cone: its facets are narrower than the bends between them".into(),
            ));
        }
    }
    let (mut outline, mut arcs) = (Vec::new(), Vec::new());
    for j in 0..n {
        let [lb, rb, ..] = corners(j);
        outline.push(lb);
        arcs.push(None);
        if let Some([bottom, _]) = cuts.get(j) {
            outline.push(rb);
            arcs.push(end(*bottom));
        }
    }
    outline.push(corners(n - 1)[1]);
    arcs.push(None);
    for j in (0..n).rev() {
        let [.., rt, lt] = corners(j);
        outline.push(rt);
        arcs.push(None);
        outline.push(lt);
        arcs.push(if j > 0 { end(cuts[j - 1][1]) } else { None });
    }
    for j in 0..n {
        layout.cells.push(Cell {
            name: format!("facet{}", j + 1),
            polygon: corners(j).to_vec(),
            arcs: vec![None; 4],
            holes: Vec::new(),
            inner: Vec::new(),
            zone: j.checked_sub(1).map_or(Zone::Base, Zone::Flange),
        });
    }
    for (j, [bottom, top]) in cuts.iter().enumerate() {
        let (s0, s1) = (bottom.s0, bottom.s0 + allowance);
        layout.cells.push(Cell {
            name: layout.folds[j].flange.clone(),
            polygon: vec![
                bottom.point(s0),
                bottom.point(s1),
                top.point(s1),
                top.point(s0),
            ],
            arcs: vec![end(*bottom), None, end(*top), None],
            holes: Vec::new(),
            inner: Vec::new(),
            zone: Zone::Bend(j),
        });
        let mid = (s0 + s1) / 2.0;
        layout.flat.bends.push(FlatBend {
            flange: layout.folds[j].flange.clone(),
            order: j + 1,
            line: [bottom.point(mid), top.point(mid)],
            angle: angle.to_degrees(),
            inner_radius: r,
            k_factor: k,
            up: false,
        });
    }
    // A cylinder's facets line up in the flat.
    crate::flat::simplify(&mut outline, &mut arcs, true);
    layout.flat.outline = outline;
    layout.flat.outline_arcs = arcs;
    Ok((layout, place))
}

/// A bend's end: straight across a cylinder's bends.
fn end(cut: Cut) -> Option<Curve> {
    (cut.c[1].abs() > EPS || cut.c[2].abs() > EPS).then_some(Curve::Cut(cut))
}

fn sub(a: P3, b: P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: P3, s: f64) -> P3 {
    a.map(|x| x * s)
}

fn dot(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: P3, b: P3) -> P3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: P3) -> f64 {
    dot(a, a).sqrt()
}

fn unit(a: P3) -> P3 {
    scale(a, 1.0 / norm(a))
}

fn add2(a: P2, b: P2) -> P2 {
    [a[0] + b[0], a[1] + b[1]]
}

fn sub2(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}

fn dot2(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn scale2(a: P2, s: f64) -> P2 {
    [a[0] * s, a[1] * s]
}

fn unit2(a: P2) -> P2 {
    scale2(a, 1.0 / a[0].hypot(a[1]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::Radial::{self, *};

    fn part(
        radii: [f64; 2],
        height: f64,
        sweep: f64,
        material: Radial,
        facets: Facets,
    ) -> PartSpec {
        let mut spec: PartSpec = serde_json::from_str(
            r#"{"id": "facets", "thickness": 2, "bend_radius": 2, "k_factor": 0.4, "base": {}}"#,
        )
        .unwrap();
        spec.base.cone = Some(Cone {
            radii,
            height,
            sweep,
            material,
            facets: Some(facets),
        });
        spec
    }

    /// Flat area, volume and surface area as a press-brake lofted flange
    /// between the same arcs measures them.
    #[test]
    fn faceted_cones_measure_as_the_reference_parts() {
        use Facets::*;
        let cases = [
            (
                part([50.0, 30.0], 60.0, 180.0, Outside, ChordTolerance(2.0)),
                6,
                7950.112048,
                15963.198117,
                16722.595362,
            ),
            (
                part([50.0, 30.0], 60.0, 180.0, Inside, ChordTolerance(2.0)),
                6,
                7628.573017,
                15320.120055,
                16059.159734,
            ),
            (
                part([30.0, 50.0], 60.0, 120.0, Outside, FacetDistance(13.0)),
                9,
                5371.424022,
                10787.548599,
                11381.881069,
            ),
            (
                part([50.0, 10.0], 20.0, 120.0, Outside, FacetAngle(25.0)),
                3,
                2643.804821,
                5299.155433,
                5726.630820,
            ),
            (
                part([50.0, 50.0], 60.0, 180.0, Outside, ChordTolerance(1.0)),
                8,
                9491.976800,
                19049.927045,
                19924.924613,
            ),
            (
                part([30.0, 50.0], 60.0, 120.0, Outside, FacetAngle(61.0)),
                1,
                4214.261501,
                8428.523002,
                8958.633344,
            ),
        ];
        for (spec, n, flat, volume, area) in cases {
            let what = format!("{:?}", spec.base.cone);
            let f = crate::flat::flat_pattern(&spec).unwrap();
            assert_eq!(f.bends.len(), n - 1, "{what}");
            assert!((f.area() - flat).abs() < 1e-5, "{what}: flat {}", f.area());
            let m = crate::fold::exact_mass(&spec).unwrap();
            assert!(
                (m.volume - volume).abs() < 1e-5,
                "{what}: volume {}",
                m.volume
            );
            assert!((m.area - area).abs() < 1e-5, "{what}: area {}", m.area);
            let mesh = crate::fold::mass(&crate::fold::fold(&spec).unwrap(), 0.01).unwrap();
            assert!(
                (mesh.volume - volume).abs() < 1e-4 * volume,
                "{what}: mesh {}",
                mesh.volume
            );
        }
    }

    /// The facet counts the reference CAD gives over sweeps of 180 (r 50 to
    /// 30) and 120 degrees (r 30 to 50, and 50 to 10).
    #[test]
    fn the_fewest_facets_within_the_limit_on_the_larger_arc() {
        use Facets::*;
        let probed = [
            (ChordTolerance(0.5), [12, 8]),
            (ChordTolerance(1.0), [8, 6]),
            (ChordTolerance(2.0), [6, 4]),
            (ChordTolerance(3.5), [5, 3]),
            (FacetAngle(10.0), [9, 6]),
            (FacetAngle(25.0), [4, 3]),
            (FacetAngle(40.0), [3, 2]),
            (FacetAngle(61.0), [2, 1]),
            (FacetDistance(5.0), [32, 21]),
            (FacetDistance(13.0), [13, 9]),
            (FacetDistance(20.0), [8, 6]),
            (FacetDistance(33.0), [5, 4]),
        ];
        for (rule, counts) in probed {
            let got = [180.0f64, 120.0].map(|s| rule.count(s.to_radians(), 50.0).unwrap());
            assert_eq!(got, counts, "{rule:?}");
        }
        assert_eq!(ChordTolerance(0.0).count(1.0, 50.0), None);
    }
}
