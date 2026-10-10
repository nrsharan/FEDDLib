#!/usr/bin/env python3
"""Writes the input files of every run of the amlodipine study into cases/.

Each run starts from the study's inputs for its artery in base/ (the files the
paper's results were computed with) and changes only what is listed here:

* every run: adaptive time stepping (see adaptive_segments), output by time,
  checkpoints, and the full stress tensor among the postprocessing fields;
* phini: the compressibility Alpha2 of the other arteries;
* narula: the remeshed geometry narula_r0.1_L1 (smoothed cap shoulders);
* each run: the parameters of its study (the table in CASES below).

cases/index.tsv lists every run with its job name, mesh and changes.
Run it from this directory; it rewrites cases/ completely.
"""

import os
import shutil
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "base")
CASES_DIR = os.path.join(HERE, "cases")

GEOMETRIES = {
    "dan": ("simulationParameters_dan.xml", "materialParameters_dan.xml", "Artery_dan_SCI.mesh"),
    "kim": ("simulationParameters_kim.xml", "materialParameters_kim_guzman.xml", "Artery_kim_guzman_SCI.mesh"),
    "narula": ("simulationParameters_narula.xml", "materialParameters_narula.xml", "narula_r0.1_L1.mesh"),
    "phini": ("simulationParameters_phini.xml", "materialParameters_phini.xml", "Artery_phinikaridou_etal-coronoary_tcfa_medium_SCI.mesh"),
    "plass": ("simulationParameters_plass.xml", "materialParameters_plasschaert_heeneman_daemem.xml", "Artery_plasschaert_heeneman_daemem_SCI.mesh"),
}

# Alpha2 by region (Volume Flag) as in dan, kim and plass; phini had 10.0 everywhere
PHINI_ALPHA2 = {15: "198.654", 16: "151.73775", 17: "151.73775", 18: "500.0",
                19: "151.73775", 20: "151.73775", 21: "151.73775", 22: "165.528"}

STRESS_FIELDS = ["Sxx", "Sxy", "Sxz", "Syx", "Syy", "Syz", "Szx", "Szy", "Szz"]

# Phase ends: load ramp, reorientation, accelerated active, reorientation, active,
# growth, reorientation, drug inflow, pressure drop, end of the run
PHASE_ENDS = [1., 20., 220., 240., 540., 840., 860., 1000., 1200., 1500.]
EXPORT_INTERVAL = 50.
CHECKPOINT_TIMES = [220., 540., 860., 1200.]
GROWTH_START = 540.
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
    segments.sort()

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
    set_parameter(simulation, stepping, "Checkpointing", "bool", "true")
    set_parameter(simulation, stepping, "Checkpoint directory", "string", "checkpoints")
    set_parameter(simulation, stepping, "Checkpoint Times", "Array(double)", "{%s}" % ", ".join(fmt(t) for t in CHECKPOINT_TIMES if t < final_time))
    set_parameter(simulation, stepping, "Export history", "bool", "true")


def common(simulation, final_time, extra_segments=(), extra_export_times=()):
    root = simulation.getroot()
    adaptive_segments(simulation, final_time, extra_segments)

    fields = parameter(sublist(root, "Parameter"), "Post Processing Fields")
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

# Revision 3.4: remeshed narula and phini geometries (original geometry refined, smoothed shoulders)
for artery in ["narula", "phini"]:
    for variant in ["original_L1", "original_L2", "r0.1_L2", "r0.2_L2"]:
        name = "%s_%s.mesh" % (artery, variant)
        add("mesh_variation", "artery_%s_%s" % (artery, variant), artery, chain(mesh(name), drug(False)), "control, mesh " + name)
        add("mesh_variation", "artery_%s_%s_withDrug" % (artery, variant), artery, mesh(name), "drug, mesh " + name)

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
                       "dose_variation": "dose", "drug_response_variation": "resp", "mesh_variation": "mesh",
                       "axial_stretch": "ax", "basal_tone_tuning": "eta"}


def job_name(study, case):
    return STUDY_ABBREVIATIONS[study] + "_" + case.replace("artery_", "").replace("withDrug", "d")


def main():
    if os.path.isdir(CASES_DIR):
        shutil.rmtree(CASES_DIR)
    rows = ["case\tjob\tartery\tmesh\tfinal time\tchanges"]
    names = set()
    for study, case, artery, change, description in CASES:
        sim_file, mat_file, mesh_file = GEOMETRIES[artery]
        simulation = read(os.path.join(BASE, sim_file))
        material_tree = read(os.path.join(BASE, mat_file))

        parameter(sublist(simulation.getroot(), "Mesh Partitioner"), "Mesh 1 Name").set("value", mesh_file)
        if artery == "phini":
            for flag, region in regions(material_tree):
                parameter(region, "Alpha2").set("value", PHINI_ALPHA2[flag])
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

        folder = os.path.join(CASES_DIR, study, case)
        os.makedirs(folder)
        write(simulation, os.path.join(folder, "simulationParameters.xml"), os.path.join(BASE, sim_file))
        write(material_tree, os.path.join(folder, mat_file), os.path.join(BASE, mat_file))

        name = job_name(study, case)
        assert name not in names, name
        names.add(name)
        mesh_name = parameter(sublist(simulation.getroot(), "Mesh Partitioner"), "Mesh 1 Name").get("value")
        with open(os.path.join(folder, "case.env"), "w") as f:
            f.write("JOB_NAME=%s\nMATERIAL=%s\nMESH=%s\n" % (name, mat_file, mesh_name))
        rows.append("\t".join([study + "/" + case, name, artery, mesh_name, "%g" % final_time, description]))
    with open(os.path.join(CASES_DIR, "index.tsv"), "w") as f:
        f.write("\n".join(rows) + "\n")
    print("%d runs written to %s" % (len(CASES), CASES_DIR))


if __name__ == "__main__":
    main()
