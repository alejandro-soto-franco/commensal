//! PyO3 bindings for `commensal`.
//!
//! One flat extension module, `_core`, that the `commensal` Python package
//! wraps in one submodule per Rust module. A struct crosses as a dict rather
//! than a Python class, since every quantity here is a handful of named
//! floats and a dict keeps the boundary thin.

use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;

#[allow(dead_code)]
fn to_py(e: commensal::Error) -> PyErr {
    match e {
        commensal::Error::Io(_) => PyOSError::new_err(e.to_string()),
        _ => PyValueError::new_err(e.to_string()),
    }
}

#[pyfunction]
fn amplicon_dereplicate(reads: Vec<Vec<u8>>) -> Vec<(Vec<u8>, u32)> {
    commensal::amplicon::dereplicate(reads)
        .into_iter()
        .map(|u| (u.sequence, u.abundance))
        .collect()
}

#[pyfunction]
#[pyo3(signature = (uniques, min_fold=2.0, min_parent=8, max_shift=16))]
fn amplicon_flag_bimeras(
    uniques: Vec<(Vec<u8>, u32)>,
    min_fold: f64,
    min_parent: u32,
    max_shift: isize,
) -> Vec<bool> {
    let uniques: Vec<commensal::Unique> = uniques
        .into_iter()
        .map(|(sequence, abundance)| commensal::Unique { sequence, abundance })
        .collect();
    let options = commensal::amplicon::BimeraOptions {
        min_fold_parent_over_abundance: min_fold,
        min_parent_abundance: min_parent,
        max_shift,
    };
    commensal::amplicon::flag_bimeras(&uniques, &options)
}

/// A struct crosses as a dict, so Python sees named fields rather than a tuple
/// whose order a caller has to remember.
#[pyfunction]
fn consortium_steady_state(py: Python<'_>, input: f64) -> PyResult<Py<PyAny>> {
    let s = commensal::consortium::Consortium::default().steady_state(input);
    let d = pyo3::types::PyDict::new(py);
    d.set_item("substrate", s.substrate)?;
    d.set_item("biomass", s.biomass)?;
    d.set_item("butyrate", s.butyrate)?;
    d.set_item("persists", s.persists)?;
    Ok(d.into_any().unbind())
}

/// A trajectory as two numpy arrays: times, and a three-column state.
type Trajectory<'py> = (Bound<'py, numpy::PyArray1<f64>>, Bound<'py, numpy::PyArray2<f64>>);

/// A trajectory crosses as two numpy arrays, times and a three-column state,
/// since a list of tuples that long is slow to build and slow to use.
#[pyfunction]
fn consortium_simulate<'py>(
    py: Python<'py>,
    dose: f64,
    hours: f64,
    dt: f64,
) -> PyResult<Trajectory<'py>> {
    use numpy::{IntoPyArray, ndarray::Array2};
    let hours_positive = hours > 0.0;
    if !hours_positive {
        return Err(PyValueError::new_err("hours must be positive"));
    }
    let dt_positive = dt > 0.0;
    if !dt_positive {
        return Err(PyValueError::new_err("dt must be positive"));
    }
    let c = commensal::consortium::Consortium::default();
    let initial = commensal::consortium::State {
        substrate: 0.0,
        biomass: 0.5,
        butyrate: 0.0,
    };
    let trace = commensal::consortium::simulate(&c, initial, move |_| dose, hours, dt);
    let times: Vec<f64> = trace.iter().map(|(t, _)| *t).collect();
    let mut states = Array2::zeros((trace.len(), 3));
    for (i, (_, s)) in trace.iter().enumerate() {
        states[[i, 0]] = s.substrate;
        states[[i, 1]] = s.biomass;
        states[[i, 2]] = s.butyrate;
    }
    Ok((times.into_pyarray(py), states.into_pyarray(py)))
}

#[pyfunction]
fn receptor_decompose(py: Python<'_>, observed_fold: f64) -> PyResult<Py<PyAny>> {
    let constructs = commensal::receptor::measured_constructs();
    let stim = constructs
        .iter()
        .find(|c| c.stimulatory)
        .ok_or_else(|| PyValueError::new_err("no stimulatory construct"))?;
    let silent = constructs
        .iter()
        .find(|c| !c.stimulatory)
        .ok_or_else(|| PyValueError::new_err("no silent construct"))?;
    let d = commensal::receptor::decompose(observed_fold, stim, silent);
    let out = pyo3::types::PyDict::new(py);
    out.set_item("observed_fold", d.observed_fold)?;
    out.set_item("monomer_fold", d.monomer_fold)?;
    out.set_item("dimer_fold", d.dimer_fold)?;
    out.set_item("proofreading_steps", d.proofreading_steps)?;
    Ok(out.into_any().unbind())
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(amplicon_dereplicate, m)?)?;
    m.add_function(wrap_pyfunction!(amplicon_flag_bimeras, m)?)?;
    m.add_function(wrap_pyfunction!(consortium_steady_state, m)?)?;
    m.add_function(wrap_pyfunction!(consortium_simulate, m)?)?;
    m.add_function(wrap_pyfunction!(receptor_decompose, m)?)?;
    Ok(())
}
