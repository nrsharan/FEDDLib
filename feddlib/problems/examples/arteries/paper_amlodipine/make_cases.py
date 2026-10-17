#!/usr/bin/env python3
"""Writes the input files of every run of the amlodipine study into cases/.

Each run starts from the study's inputs for its artery in base/ (the files the
paper's results were computed with) and changes only what is listed here:

* every run: adaptive time stepping (see adaptive_segments), output by time and
  checkpoints;
* phini: Alpha2 of 100.0 (PHINI_ALPHA2), an adaptive load ramp that follows the time
  step (with fixed load steps its Newton iteration diverges, see PHINI_RAMP_DT), smaller minimum
  time steps (PHINI_MIN_DT_DIVISOR) and more linear iterations (PHINI_GMRES_ITERATIONS, the case's
  own solverParameters.xml);
* every run: the stress tensor (Sxx ... Szz) and the fibre stresses (SFiber1, SFiber2)
  among the postprocessing fields;
* every run: the preconditioner is built in the first PRECONDITIONER_NEWTON_STEPS Newton
  steps of each time step and reused in the later ones;
* every run: growthDirectionBool GROWTH_DIRECTION_BOOL in every region (growth directions
  kept at their values at the start of growth);
* kim, narula, phini: the paper's geometry remeshed finer at the cap (KIM_MESH, NARULA_MESH,
  PHINI_MESH); on hold until the mesh convergence study has run on the corrected meshes;
* each run: the parameters of its study (the table in CASES below).

cases/index.tsv lists every run with its job name, mesh, nodes and changes. The runs of
kim, narula and phini are on hold (HOLD in case.env; submit_cases.sh skips them): their mesh is
decided by the mesh convergence study, which has to be repeated on the corrected meshes.
Run it from this directory; it rewrites cases/ completely.
"""

import os
import shutil
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "base")
CASES_DIR = os.path.join(HERE, "cases")
MESHES = os.path.join(HERE, "..", "..", "..", "..", "..", "meshes", "SPP2311")

# kim: the paper's geometry (Artery_kim_guzman_SCI.mesh, cap elements ~0.088 mm) remeshed with the
# pipeline of the mesh convergence study (cross-section recovered from the paper's mesh, Gmsh, global
# size 0.09 mm, 0.03 mm near the cap shoulders and stress peaks; cap elements ~0.054 mm)
KIM_MESH = "kim_original_L2.mesh"
# narula, phini: level L2 of the mesh convergence study, unsmoothed shoulders (the paper's geometry).
# The remeshed meshes written before 2026-10-16 had half of their boundary triangles pointing into
# the volume; FEDDLib takes a surface load's normal from the vertex order, so the lumen pressure
# cancelled and every run on them (the first mesh convergence study, kim L2) was unloaded. The
# meshes are now written with outward triangles; KIM_MESH, NARULA_MESH and PHINI_MESH are to be
# decided again from the repeated study (HELD)
NARULA_MESH = "narula_original_L2.mesh"
PHINI_MESH = "phini_original_L2.mesh"

GEOMETRIES = {
    "dan": ("simulationParameters_dan.xml", "materialParameters_dan.xml", "Artery_dan_SCI.mesh"),
    "kim": ("simulationParameters_kim.xml", "materialParameters_kim_guzman.xml", KIM_MESH),
    "narula": ("simulationParameters_narula.xml", "materialParameters_narula.xml", NARULA_MESH),
    "phini": ("simulationParameters_phini.xml", "materialParameters_phini.xml", PHINI_MESH),
    "plass": ("simulationParameters_plass.xml", "materialParameters_plasschaert_heeneman_daemem.xml", "Artery_plasschaert_heeneman_daemem_SCI.mesh"),
}

# phini's Alpha2 in every region: 100.0. With the values of dan, kim and plass (198.654, 151.73775,
# 151.73775, 500.0, 151.73775, 151.73775, 151.73775, 165.528 in regions 15-22) phini original stopped on
# the corrected remeshes: L1 at 0.91 s and L2 at 52.4 s, its Newton iteration stalling at any time step
# (mesh convergence study, 2026-10-17). With the paper's 10.0 (FEDDLib-Sharan 18a1448b) the paper's cap
# lost ~10 % of its volume elastically (DetFe p50 0.88-0.90, nodal min 0.72; dan's cap: 0.99). phini
# original L2 ran to 1500 s with 50.0 (cap DetFe p50 0.957) and with 100.0 (2026-10-17)
PHINI_ALPHA2 = "100.0"
# phini's segments after the load ramp may shrink to their dt / PHINI_MIN_DT_DIVISOR (the other
# arteries: to their dt)
PHINI_MIN_DT_DIVISOR = 16
# phini's GMRES: at most PHINI_GMRES_ITERATIONS iterations without restart (solverParameters.xml: 1000).
# With 1000 it missed its tolerance in most Newton steps of the phini runs
PHINI_GMRES_ITERATIONS = "3000"

# phini's load ramp adapts its time step. With the consistent Alpha2 and fixed load steps its first
# Newton iteration diverges at 0.02 s on the paper's mesh and at 0.01 s on the remeshes; at 0.005 s the
# finer remeshes (original L1, L2, r0.1 L2) still diverge at 84-92 % of the pressure. The ramp follows
# the time step (main.cpp, "Ramp Follows Time Step"), starts with PHINI_RAMP_DT and may shrink to
# PHINI_RAMP_MIN_DT (phini original L2 needed 0.00039 at 84 %, 2026-10-17). FEDDLib requires
# Load Step Size = the first segment's dt.
PHINI_RAMP_DT = "0.01"
PHINI_RAMP_MIN_DT = "0.0001"
PHINI_RAMP_MAX_DT = "0.05"

STRESS_FIELDS = ["Sxx", "Sxy", "Sxz", "Syx", "Syy", "Syz", "Szx", "Szy", "Szz", "SFiber1", "SFiber2"]

# Remeshed meshes: tetrahedra per node (48 ranks)
TETS_PER_NODE = {"narula": 25000, "phini": 30000, "kim": 30000}
# Meshes whose node count is set by a benchmark instead. narula L3 (1.85M unknowns) on 13 nodes needed
# ~20 s per Newton step, twice L2's at the same load per rank: FROSch's coarse problem grows with the
# number of ranks and is solved on one. 3 load steps of narula original L3 (2026-10-17,
# runs_frozen/bench_narula_L3): 13 nodes ~20 s per Newton step (0.072 node-h), 7 nodes 22.4 s (0.044),
# 4 nodes 31 s (0.034)
MESH_NODES = {"narula_original_L3.mesh": 4, "narula_r0.1_L3.mesh": 4}
# The L3 meshes (5-12 MB) are not in the repository (smoothing_package/<artery>/meshes on the
# user's machine); their sizes, for the node counts when the file is not there
MESH_TETRAHEDRA = {"narula_original_L3.mesh": 316461, "narula_r0.1_L3.mesh": 316935,
                   "phini_original_L3.mesh": 147287, "phini_r0.1_L3.mesh": 147911}
HELD = "the repeated mesh convergence study decides the mesh of kim, narula and phini"

# Phase ends: load ramp, reorientation, accelerated active, reorientation, active,
# growth, reorientation, drug inflow, pressure drop, end of the run
PHASE_ENDS = [1., 20., 220., 240., 540., 840., 860., 1000., 1200., 1500.]
EXPORT_INTERVAL = 50.
# Successful time steps of an adaptive segment before its time step grows (FEDDLib's default: 5)
CONVERGED_STEPS_BEFORE_INCREASE = 7
GROWTH_START = 540.
GROWTH_END = 840.
# Checkpoints at the phase ends and every 25 s of growth, so a run can be repeated from
# shortly before the time of interest
GROWTH_CHECKPOINT_INTERVAL = 25.
CHECKPOINT_TIMES = sorted([220., 540., 860., 1200.] + [GROWTH_START + GROWTH_CHECKPOINT_INTERVAL * k for k in range(1, int((GROWTH_END - GROWTH_START) / GROWTH_CHECKPOINT_INTERVAL) + 1)])
# Newton steps of a time step that build the preconditioner; the later ones reuse it
# (the paper: rebuilt in every Newton step)
PRECONDITIONER_NEWTON_STEPS = 2
# growthDirectionBool of every region: 0 keeps the growth directions at their values at the start of
# growth; with 1 (the paper) the element takes them from the Mandel stress in every step, and in the
# fibrous cap they flip between two orientations from one step to the next (a residual at the start of
# each step that does not shrink with dt; narula r0.1_L1 and plass could not pass ~595 and ~721 s)
GROWTH_DIRECTION_BOOL = "0"
NO_DRUG = "10000000.0"  # Inflow Start Time of the drug-free runs (past the end)


def fmt(x):
    return repr(float(x))


def read(path):
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    return ET.parse(path, parser)


def write(tree, path, base):
    """Writes the file in the style of its input (<Parameter .../> or <Parameter ... />)."""
    text = ET.tostring(tree.getroot(), encoding="unicode")
    if '"/>' in open(base).read():
        text = text.replace('" />', '"/>')
    with open(path, "w") as f:
        f.write(text + "\n")


def sublist(element, *names):
    for name in names:
        found = element.find("ParameterList[@name='%s']" % name)
        assert found is not None, (name, names)
        element = found
    return element


def parameter(element, name):
    found = element.find("Parameter[@name='%s']" % name)
    assert found is not None, name
    return found


def indent(depth):
    return "\n" + "    " * depth


def depth_of(tree, element):
    parents = {c: p for p in tree.iter() for c in p}
    depth = 0
    while element in parents:
        element = parents[element]
        depth += 1
    return depth


def set_parameter(tree, element, name, ptype, value):
    """Sets a parameter, or appends it after the last parameter of the list."""
    found = element.find("Parameter[@name='%s']" % name)
    if found is None:
        found = ET.Element("Parameter", {"name": name, "type": ptype, "value": value})
        children = list(element)
        last = max(i for i, c in enumerate(children) if c.tag == "Parameter")
        while last + 1 < len(children) and children[last + 1].tag is ET.Comment and "\n" not in (children[last].tail or ""):
            last += 1  # a comment on the same line belongs to the parameter
        found.tail = children[last].tail
        children[last].tail = indent(depth_of(tree, element) + 1)
        element.insert(last + 1, found)
    found.set("type", ptype)
    found.set("value", value)


def regions(material):
    for region in sublist(material.getroot(), "Parameter Solid").findall("ParameterList"):
        yield int(parameter(region, "Volume Flag").get("value")), region


def set_material(material, name, value, flags=None):
    for flag, region in regions(material):
        if flags is None or flag in flags:
            parameter(region, name).set("value", value)


def scale_material(material, name, factor):
    for _, region in regions(material):
        p = parameter(region, name)
        p.set("value", fmt(float(p.get("value")) * factor))


def adaptive_segments(simulation, final_time, extra_segments=()):
    """Every segment but the load ramp adapts its time step: it starts with the paper's dt, which is
    also its smallest one, and may grow to half the segment's length. Elements that do not converge
    locally are accepted at the smallest time step in the growth segment only."""
    stepping = sublist(simulation.getroot(), "Timestepping Parameter")
    intervals = sublist(stepping, "Timestepping Intervalls")
    segments = []
    for item in intervals.findall("ParameterList"):
        segments.append((float(parameter(item, "Start Time").get("value")), parameter(item, "dt").get("value")))
    segments += list(extra_segments)
    # a run that ends early (basal_tone_tuning) leaves out the segments it does not reach
    segments = sorted(s for s in segments if s[0] < final_time)

    depth = depth_of(simulation, intervals)
    for child in list(intervals):
        if child.tag == "ParameterList":
            intervals.remove(child)
    number = parameter(intervals, "Number of Segments")
    number.set("value", str(len(segments)))
    number.tail = indent(depth + 1)
    for i, (start, dt) in enumerate(segments):
        end = segments[i + 1][0] if i + 1 < len(segments) else final_time
        item = ET.SubElement(intervals, "ParameterList", {"name": str(i + 1)})
        item.text = indent(depth + 2)
        item.tail = indent(depth + 1)
        values = [("Start Time", "double", fmt(start)), ("dt", "double", dt)]
        if i == 0:
            values.append(("Adaptive", "bool", "false"))
        else:
            values += [("Adaptive", "bool", "true"),
                       ("Minimum dt", "double", dt),
                       ("Maximum dt", "double", fmt(0.5 * (end - start))),
                       ("Accept Element Failures At Minimum dt", "bool", "true" if start == GROWTH_START else "false")]
        for name, ptype, value in values:
            p = ET.SubElement(item, "Parameter", {"name": name, "type": ptype, "value": value})
            p.tail = indent(depth + 2)
        p.tail = indent(depth + 1)
    item.tail = indent(depth)

    set_parameter(simulation, stepping, "Final time", "double", fmt(final_time))
    set_parameter(simulation, stepping, "Adaptive Time Stepping", "bool", "false")
    set_parameter(simulation, stepping, "Converged Time Steps Before Increase", "int", str(CONVERGED_STEPS_BEFORE_INCREASE))
    set_parameter(simulation, stepping, "Checkpointing", "bool", "true")
    set_parameter(simulation, stepping, "Checkpoint directory", "string", "checkpoints")
    set_parameter(simulation, stepping, "Checkpoint Times", "Array(double)", "{%s}" % ", ".join(fmt(t) for t in CHECKPOINT_TIMES if t < final_time))
    set_parameter(simulation, stepping, "Export history", "bool", "true")


def common(simulation, final_time, extra_segments=(), extra_export_times=()):
    root = simulation.getroot()
    adaptive_segments(simulation, final_time, extra_segments)

    p = sublist(root, "Parameter")
    set_parameter(simulation, p, "Rebuild Preconditioner every Newton Iteration", "bool", "false")
    set_parameter(simulation, p, "newtonLimit", "int", str(PRECONDITIONER_NEWTON_STEPS))

    fields = parameter(p, "Post Processing Fields")
    names = [f.strip() for f in fields.get("value").strip("{}").split(",")]
    fields.set("value", "{%s}" % ", ".join(names + [f for f in STRESS_FIELDS if f not in names]))

    exporter = sublist(root, "Exporter")
    times = sorted(set([t for t in PHASE_ENDS if t <= final_time] + list(extra_export_times)))
    set_parameter(simulation, exporter, "Export Times", "Array(double)", "{%s}" % ", ".join(fmt(t) for t in times))
    set_parameter(simulation, exporter, "Export Interval", "double", fmt(EXPORT_INTERVAL))


# ---------------------------------------------------------------------------------------------
# The runs. Each entry: (study folder, case folder, artery, change function, description)
# A change function gets (simulation, material) and returns the final time and, optionally,
# extra time segments and output times.

def drug(on):
    def change(sim, mat):
        if not on:
            parameter(sublist(sim.getroot(), "Parameter"), "Inflow Start Time").set("value", NO_DRUG)
    return change


def chain(*changes):
    def change(sim, mat):
        for c in changes:
            c(sim, mat)
    return change


def unloading(sim, mat):
    """Pressure from 85 mmHg to 0 over 1500-1600 s, then held to 1700 s (residual stresses)."""
    p = sublist(sim.getroot(), "Parameter")
    parameter(p, "Pressure Reduction Start Time").set("value", "1500.0")
    parameter(p, "Pressure Reduction End Time").set("value", "1600.0")
    parameter(p, "Pressure Reduction Amount mmHg").set("value", "85.0")


def pressure_drop(mmhg):
    def change(sim, mat):
        parameter(sublist(sim.getroot(), "Parameter"), "Pressure Reduction Amount mmHg").set("value", fmt(mmhg))
    return change


def material(name, value, flags=None):
    return lambda sim, mat: set_material(mat, name, value, flags)


def scaled(name, factor):
    return lambda sim, mat: scale_material(mat, name, factor)


def inflow_concentration(c):
    return lambda sim, mat: set_parameter(sim, sublist(sim.getroot(), "Parameter"), "Inflow Concentration", "double", fmt(c))


def axial_stretch(s):
    return lambda sim, mat: set_parameter(sim, sublist(sim.getroot(), "Parameter"), "Axial Stretch", "double", fmt(s))


def phini_adaptive_ramp(simulation):
    """phini's load ramp (the first time segment) adapts its time step and the load follows it."""
    set_parameter(simulation, sublist(simulation.getroot(), "Parameter"), "Ramp Follows Time Step", "bool", "true")
    ramp = sublist(simulation.getroot(), "Timestepping Parameter", "Timestepping Intervalls", "1")
    assert float(parameter(ramp, "Start Time").get("value")) == 0. and parameter(ramp, "dt").get("value") == PHINI_RAMP_DT
    set_parameter(simulation, ramp, "Adaptive", "bool", "true")
    set_parameter(simulation, ramp, "Minimum dt", "double", PHINI_RAMP_MIN_DT)
    set_parameter(simulation, ramp, "Maximum dt", "double", PHINI_RAMP_MAX_DT)


def phini_minimum_dt(simulation):
    """phini's segments after the load ramp may shrink to their dt / PHINI_MIN_DT_DIVISOR."""
    for item in sublist(simulation.getroot(), "Timestepping Parameter", "Timestepping Intervalls").findall("ParameterList"):
        if float(parameter(item, "Start Time").get("value")) > 0.:
            parameter(item, "Minimum dt").set("value", fmt(float(parameter(item, "dt").get("value")) / PHINI_MIN_DT_DIVISOR))


def phini_solver(path):
    """solverParameters.xml with PHINI_GMRES_ITERATIONS iterations and Krylov vectors."""
    solver = read(os.path.join(HERE, "solverParameters.xml"))
    gmres = sublist(solver.getroot(), "ThyraSolver", "Linear Solver Types", "Belos", "Solver Types", "Block GMRES")
    parameter(gmres, "Maximum Iterations").set("value", PHINI_GMRES_ITERATIONS)
    parameter(gmres, "Num Blocks").set("value", PHINI_GMRES_ITERATIONS)
    write(solver, path, os.path.join(HERE, "solverParameters.xml"))


def mesh(name):
    return lambda sim, mat: parameter(sublist(sim.getroot(), "Mesh Partitioner"), "Mesh 1 Name").set("value", name)


CASES = []
UNLOADING = {}  # case -> True: runs to 1700 s with the unloading


def add(study, case, artery, change, description, unload=False):
    CASES.append((study, case, artery, change, description))
    if unload:
        UNLOADING[(study, case)] = True


# Paper: kappa (nitric oxide) variation, dan; the 30 set is the baseline
for level, media, degen in [("15", "126.023", "84.4352"), ("30", None, None), ("45", "81.5441", "54.6345")]:
    kappa = (lambda m, d: chain(material("Kappa", m, {16}), material("Kappa", d, {17})))(media, degen) if media else drug(True)
    add("kappa_variation", "artery_dan_" + level, "dan", chain(kappa, drug(False)),
        "control, Kappa media/degenerated media %s" % ("103.783/69.5349" if media is None else media + "/" + degen),
        unload=(level == "30"))
    add("kappa_variation", "artery_dan_withDrug_" + level, "dan", kappa,
        "drug, Kappa media/degenerated media %s" % ("103.783/69.5349" if media is None else media + "/" + degen))

# Paper: the other four arteries, control and drug (dan: kappa_variation 30)
for artery in ["kim", "narula", "phini", "plass"]:
    add("arteries", "artery_" + artery, artery, drug(False), "control", unload=True)
    add("arteries", "artery_%s_withDrug" % artery, artery, drug(True), "drug")

# Paper: pressure drop of 10, 20, 30 mmHg over 1000-1200 s, drug (0 mmHg: kappa_variation 30 and arteries)
for mmhg in [10, 20, 30]:
    for artery in ["dan", "kim", "narula", "phini", "plass"]:
        add("pressure_drop_variation/%d_pressure_drop" % mmhg, "artery_%s_withDrug" % artery, artery,
            pressure_drop(mmhg), "drug, pressure drop %d mmHg" % mmhg)

# Paper: diffusivity D0 of adventitia (15) and degenerated media (17); e = baseline
for letter, adventitia, degen in [("a", "0.0007", "0.00035"), ("b", "0.0007", "0.0035"), ("c", "0.0007", "0.035"),
                                  ("d", "0.007", "0.00035"), ("f", "0.007", "0.035"),
                                  ("h", "0.07", "0.00035"), ("i", "0.07", "0.0035"), ("j", "0.07", "0.035")]:
    add("diffusion_variation", "artery_dan_withDrug_" + letter, "dan",
        chain(material("D0", adventitia, {15}), material("D0", degen, {17})),
        "drug, D0 adventitia/degenerated media %s/%s" % (adventitia, degen))

# Paper: reaction M of media (16) and degenerated media (17); the sign follows the element in use
for letter, media, degen in [("a", "-3.e-1", "-1.5e-1"), ("b", "-3.e-2", "-1.5e-2"), ("c", "-3.e-3", "-1.5e-3"),
                             ("d", "-3.e-4", "-1.5e-4"), ("e", "-3.e-5", "-1.5e-5"), ("f", "-3.e-6", "-1.5e-6")]:
    add("reaction_variation", "artery_dan_withDrug_" + letter, "dan",
        chain(material("M", media, {16}), material("M", degen, {17})),
        "drug, M media/degenerated media %s/%s" % (media, degen))

# Revision 3.1: contractility (Kappa) of the degenerated media (17), as a share of the media's 103.783
for share, kappa in [("33", "34.25"), ("10", "10.38"), ("0", "0.0")]:
    add("degen_kappa_variation", "artery_dan_%s" % share, "dan", chain(material("Kappa", kappa, {17}), drug(False)),
        "control, Kappa degenerated media %s (%s %% of the media)" % (kappa, share))
    add("degen_kappa_variation", "artery_dan_withDrug_%s" % share, "dan", material("Kappa", kappa, {17}),
        "drug, Kappa degenerated media %s (%s %% of the media)" % (kappa, share))

# Revision 3.2: drug concentration at the walls (paper: 2 uM)
for c in [0.06, 0.2, 0.5, 1.0, 1.5]:
    add("dose_variation", "artery_dan_withDrug_%guM" % c, "dan", inflow_concentration(c), "drug, %g uM at the walls" % c)

# Revision 3.2: drug-calcium response C50 and P (all regions), 2 uM
for name, factor in [("C50", 0.5), ("C50", 2.0), ("P", 0.5), ("P", 2.0)]:
    add("drug_response_variation", "artery_dan_withDrug_%s_x%g" % (name, factor), "dan", scaled(name, factor),
        "drug, %s x%g (all regions)" % (name, factor))

# Revision 3.4: mesh convergence, narula and phini remeshed with Gmsh: the paper's geometry
# ("original") and its cap shoulders smoothed with a 0.1 mm fillet ("r0.1"); the size near the
# shoulders halves from level to level (L0 0.12, L1 0.06, L2 0.03, L3 0.015 mm)
for artery in ["phini", "narula"]:
    for family in ["original", "r0.1"]:
        for level in ["L0", "L1", "L2", "L3"]:
            name = "%s_%s_%s.mesh" % (artery, family, level)
            add("mesh_convergence", "artery_%s_%s_%s" % (artery, family, level), artery, chain(mesh(name), drug(False)), "control, mesh " + name)
            add("mesh_convergence", "artery_%s_%s_%s_withDrug" % (artery, family, level), artery, mesh(name), "drug, mesh " + name)

# Revision 3.6: 5 % axial pre-stretch, ramped with the pressure
add("axial_stretch", "artery_dan_stretch5", "dan", chain(axial_stretch(0.05), drug(False)), "control, 5 % axial stretch")
add("axial_stretch", "artery_dan_withDrug_stretch5", "dan", axial_stretch(0.05), "drug, 5 % axial stretch")

# Revision 3.3 (tuning): the active phase 240-540 s with lower Eta, to find the values that give a
# pre-drug n_C of about 0.5 and 0.3; the basal-tone runs follow once they are known
for factor in [0.75, 0.5, 0.35, 0.25, 0.15, 0.1]:
    add("basal_tone_tuning", "artery_dan_eta_x%g" % factor, "dan", chain(scaled("Eta", factor), drug(False)),
        "control to 540 s, Eta x%g (all regions)" % factor)

STUDY_ABBREVIATIONS = {"kappa_variation": "kap", "arteries": "art", "pressure_drop_variation/10_pressure_drop": "pd10",
                       "pressure_drop_variation/20_pressure_drop": "pd20", "pressure_drop_variation/30_pressure_drop": "pd30",
                       "diffusion_variation": "dif", "reaction_variation": "rea", "degen_kappa_variation": "dkap",
                       "dose_variation": "dose", "drug_response_variation": "resp", "mesh_convergence": "mesh",
                       "axial_stretch": "ax", "basal_tone_tuning": "eta"}


def job_name(study, case):
    return STUDY_ABBREVIATIONS[study] + "_" + case.replace("artery_", "").replace("withDrug", "d")


def tetrahedra(mesh_name):
    if not os.path.exists(os.path.join(MESHES, mesh_name)) and mesh_name in MESH_TETRAHEDRA:
        return MESH_TETRAHEDRA[mesh_name]
    with open(os.path.join(MESHES, mesh_name)) as f:
        words = iter(f.read().split())
    for word in words:
        if word == "Tetrahedra":
            return int(next(words))
    raise ValueError("no tetrahedra in " + mesh_name)


def main():
    if os.path.isdir(CASES_DIR):
        shutil.rmtree(CASES_DIR)
    rows = ["case\tjob\tartery\tmesh\tnodes\tfinal time\thold\tchanges"]
    names = set()
    for study, case, artery, change, description in CASES:
        sim_file, mat_file, mesh_file = GEOMETRIES[artery]
        simulation = read(os.path.join(BASE, sim_file))
        material_tree = read(os.path.join(BASE, mat_file))

        parameter(sublist(simulation.getroot(), "Mesh Partitioner"), "Mesh 1 Name").set("value", mesh_file)
        if artery == "phini":
            set_material(material_tree, "Alpha2", PHINI_ALPHA2)
            parameter(sublist(simulation.getroot(), "Parameter"), "Load Step Size").set("value", PHINI_RAMP_DT)
            ramp = sublist(simulation.getroot(), "Timestepping Parameter", "Timestepping Intervalls", "1")
            assert float(parameter(ramp, "Start Time").get("value")) == 0.
            parameter(ramp, "dt").set("value", PHINI_RAMP_DT)
        change(simulation, material_tree)

        final_time, extra_segments, extra_exports = 1500., (), ()
        if (study, case) in UNLOADING:
            unloading(simulation, material_tree)
            last_dt = sorted((float(parameter(i, "Start Time").get("value")), parameter(i, "dt").get("value"))
                             for i in sublist(simulation.getroot(), "Timestepping Parameter", "Timestepping Intervalls").findall("ParameterList"))[-1][1]
            final_time, extra_segments, extra_exports = 1700., [(1500., last_dt)], (1600., 1700.)
            description += ", unloading 1500-1600 s to 1700 s"
        if study == "basal_tone_tuning":
            final_time = 540.
        common(simulation, final_time, extra_segments, extra_exports)
        if artery == "phini":
            phini_adaptive_ramp(simulation)
            phini_minimum_dt(simulation)
        set_material(material_tree, "growthDirectionBool", GROWTH_DIRECTION_BOOL)

        folder = os.path.join(CASES_DIR, study, case)
        os.makedirs(folder)
        write(simulation, os.path.join(folder, "simulationParameters.xml"), os.path.join(BASE, sim_file))
        write(material_tree, os.path.join(folder, mat_file), os.path.join(BASE, mat_file))
        if artery == "phini":
            phini_solver(os.path.join(folder, "solverParameters.xml"))

        name = job_name(study, case)
        assert name not in names, name
        names.add(name)
        mesh_name = parameter(sublist(simulation.getroot(), "Mesh Partitioner"), "Mesh 1 Name").get("value")
        nodes = 1
        if mesh_name in MESH_NODES:
            nodes = MESH_NODES[mesh_name]
        elif artery in TETS_PER_NODE:
            nodes = -(-tetrahedra(mesh_name) // TETS_PER_NODE[artery])
        hold = HELD if artery in ("kim", "narula", "phini") and study != "mesh_convergence" else ""
        with open(os.path.join(folder, "case.env"), "w") as f:
            f.write("JOB_NAME=%s\nMATERIAL=%s\nMESH=%s\nNODES=%d\n" % (name, mat_file, mesh_name, nodes))
            if hold:
                f.write("HOLD=\"%s\"\n" % hold)
        rows.append("\t".join([study + "/" + case, name, artery, mesh_name, str(nodes), "%g" % final_time, "on hold" if hold else "", description]))
    with open(os.path.join(CASES_DIR, "index.tsv"), "w") as f:
        f.write("\n".join(rows) + "\n")
    print("%d runs written to %s" % (len(CASES), CASES_DIR))


if __name__ == "__main__":
    main()
