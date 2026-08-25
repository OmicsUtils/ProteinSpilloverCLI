Spatial protein spillover CLI

protein-spillover is an end-to-end workflow for measuring and correcting cell-segmentation-associated protein spillover in Xenium and other multiplex protein images stored in SpatialData.

The current correction model is intentionally narrow and preservation-first. It is designed primarily to improve immune-lineage protein measurements without using cell-type annotations to decide what a cell should express.

The central rule is now:

Correct only protein signal that can be physically supported by a specific cell-cell interface. If the image cannot distinguish intrinsic expression from neighbor-derived signal, preserve the measurement and flag the ambiguity.

All analyzed proteins remain available in the outputs. Only markers listed in correction_channels are altered by the spillover model.

Current correction model

The current algorithm is immune-pairwise-interface-v1.

For each selected correction marker and each directly contacting cell pair, the workflow:

identifies the focal-cell pixels belonging to that specific cell-cell interface;

assigns each focal interface pixel to at most one neighboring cell;

measures signal inside the focal interface;

establishes a localization-appropriate focal-cell reference region;

checks whether the neighboring cell is a plausible physical source of the marker;

measures only the excess signal observed inside the focal cell relative to its own reference;

converts that excess into a whole-cell-equivalent contamination amount using the fraction of valid focal-cell pixels occupied by the interface;

evaluates all configured correction scenarios from the same physically bounded interface evidence;

automatically recommends only none, conservative, or medium;

preserves and flags cases where intrinsic signal cannot be distinguished from neighbor-derived signal.

Neighbor brightness is used only to decide whether a neighbor is a plausible source. It does not directly determine how much signal is subtracted from the focal cell.

What the workflow does

The CLI runs fourteen checkpointed stages:

01_roi_selection
02_cropped_arrays
03_corrected_image
04_thresholds
05_spillover_features
06_merged_features
07_contact_graph
08_geometry_density
09_neighbor_exposure
10_correction_scenarios
11_recommendations
12_roi_h5ad
13_qc_plots
14_summary

Conceptually, the pipeline is:

SpatialData image + segmentation + table
                |
                v
native-grid registration and ROI crop
                |
                v
XOA/generic intensity preprocessing + QC masks
                |
                v
original per-cell protein measurements
                |
                v
direct cell-cell contact graph
                |
                v
pair-specific interface geometry
                |
                v
marker-localization-specific focal reference
                |
                v
source plausibility + directionality tests
                |
                v
physically supported interface contamination
                |
                +-----------------------------+
                |                             |
                v                             v
all correction scenarios              pair-level audit evidence
                |
                v
preservation-first recommendation
                |
                v
ROI AnnData + tables + QC outputs

The input Zarr remains read-only.

Core assumptions

The workflow assumes:

the table has one row per segmented biological cell;

a positive integer table column directly identifies the same cell in the label raster;

the protein image and label raster share the same full-resolution native y-by-x grid;

table spatial coordinates can be related to native image pixels using one scale, optional x/y reflection, and optional global origins;

original image-derived measurements must remain available for auditing;

correction should not require trusted cell-type labels;

large segmentation failures cannot be fully repaired by intensity subtraction.

The CLI does not infer a label crosswalk from centroid matching or polygon overlap. The configured table label column is authoritative and is checked against the cropped raster.

Files used with this workflow

protein_spillover_cli.py
protein_spillover_config.json
protein_spillover_README.md
submit_protein_spillover_all_rois.sh
run_protein_spillover_all_rois_sbatch_direct.slurm
review_spillover_phenotype_ambiguity.py   optional downstream review step

Installation

Use Python 3.10 or newer in an isolated environment containing the project's SpatialData stack.

python -m pip install -e .

Verify:

protein-spillover --help
protein-spillover run --help

First inspect a new project

protein-spillover inspect \
  --zarr /data/project/spatial_data.zarr

A more detailed inspection:

protein-spillover inspect \
  --zarr /data/project/spatial_data.zarr \
  --table cell_table \
  --image morphology_focus \
  --labels cell_labels \
  --roi-col roi

Use the output to verify table/image/label names, spatial coordinates, authoritative raster labels, ROI labels, and protein-channel names.

Configuration

Generate a template with:

protein-spillover template-config \
  --output protein_spillover_config.json

Run with:

protein-spillover run \
  --config protein_spillover_config.json

CLI values override JSON values.

Analysis channels versus correction channels

These are now deliberately separate concepts.

analysis_channels

All selected channels are measured, retained, exported, and available to reports.

For example, a panel can include lineage, state, functional, tumor, and structural proteins:

"analysis_channels": [
  "CD3E",
  "CD4",
  "CD8A",
  "PD-1",
  "LAG-3",
  "CD45RA",
  "CD45RO",
  "CD16",
  "GranzymeB",
  "CD20",
  "CD138",
  "HLA-DR",
  "CD11c",
  "CD68",
  "CD163",
  "PD-L1",
  "CD45",
  "Ki-67",
  "CD31",
  "PTEN",
  "PanCK",
  "E-Cadherin"
]

correction_channels

Only these markers are spillover-corrected. All other analyzed proteins are written through every correction scenario unchanged.

The current immune-lineage-focused set is:

"correction_channels": [
  "CD45",
  "CD3E",
  "CD4",
  "CD8A",
  "CD20",
  "CD138",
  "CD16",
  "CD11c",
  "CD68",
  "CD163",
  "HLA-DR"
]

This intentionally leaves markers such as PD-1, LAG-3, CD45RA, CD45RO, GranzymeB, Ki-67, PD-L1, PanCK, E-Cadherin, CD31, and PTEN uncorrected unless they are explicitly added later.

Repeated CLI flags can override the list:

--correction-channel CD45 \
--correction-channel CD3E \
--correction-channel CD8A

Marker localization

Every correction marker requires a localization class.

Current project mapping:

"marker_localization": {
  "CD45": "membrane",
  "CD3E": "membrane",
  "CD4": "membrane",
  "CD8A": "membrane",
  "CD20": "membrane",
  "CD138": "membrane",
  "CD16": "membrane",
  "CD11c": "membrane",
  "CD68": "intracellular",
  "CD163": "membrane",
  "HLA-DR": "membrane"
}

Accepted classes are:

membrane

intracellular

nuclear

The current pipeline has whole-cell segmentation but no true nuclear mask. A selected nuclear correction marker is therefore preserved with reference_quality = nuclear_reference_unavailable rather than being corrected using an eroded cell center as a fake nucleus.

CLI override syntax:

--marker-localization CD8A=membrane \
--marker-localization CD68=intracellular

Intensity preprocessing

xenium_xoa

Use this for Xenium Onboard Analysis protein images.

XOA images are already deconvolved, autofluorescence-background-subtracted, saturation-masked, and spectrally corrected. XOA normally stores these values with a positive offset of 100.

The workflow:

subtracts the configured XOA storage offset;

retains signed offset-adjusted values for intensity summaries;

creates a clipped nonnegative signal representation for thresholds and interface evidence;

does not perform another Gaussian background subtraction.

"input_intensity_mode": "xenium_xoa",
"xenium_xoa_intensity_offset": 100.0

apply_gaussian_background_subtraction is a legacy compatibility key and is not needed in a new Xenium XOA configuration.

precorrected

Uses supplied values directly as signed analysis intensities.

generic_gaussian

Performs broad Gaussian background subtraction for generic multiplex images that have not already been background corrected.

QC masks

A channel-matched QC mask may be supplied.

For Xenium XOA saturation masks, valid pixels are normally 0 and masked pixels are 255.

"protein_qc_mask_name": "morphology_focus_qc_masks",
"qc_mask_valid_value": 0,
"xenium_require_qc_mask": false

When no QC mask is supplied in XOA mode and xenium_zero_is_invalid_without_qc_mask is true, exact stored zeros are excluded as a conservative fallback. The official QC mask remains preferable when available.

Existing per-cell feature extraction

Stages 01-08 continue to calculate the original audit/QC features, including:

whole-cell intensity;

eroded interior intensity;

internal boundary intensity;

outer-ring intensity;

outer noncell signal;

outer other-cell signal;

positive-pixel fractions;

angular boundary coverage;

boundary anisotropy;

direct contact counts;

shared-boundary geometry;

dense-small-cell geometry scores.

These measurements remain useful for QC and downstream reporting.

The dense-small-cell geometry score no longer changes the correction amount or recommendation. It is a diagnostic only.

Pairwise interface geometry

The correction model begins from the segmentation itself.

For every touching pair A-B, the workflow constructs directed focal-cell interface bands:

A-facing-B pixels
B-facing-A pixels

The inward interface-band width is controlled by:

"interface_band_pixels": 2

or:

--interface-band 2

Each focal-cell interface pixel is assigned to one neighboring label only. This prevents one suspicious pixel from being counted multiple times when a cell contacts several neighbors.

The interface geometry is built once per ROI and reused across correction markers.

Source plausibility

For a directed pair source -> focal, the source must be genuinely marker-positive at its side of the shared interface.

A source is considered plausible only when it has enough valid interface pixels and either:

a sufficient positive-pixel fraction at that interface; or

a source-interface mean at or above the channel threshold.

The main settings are:

"minimum_interface_valid_pixels": 2,
"interface_source_positive_fraction": 0.25

The source intensity is a gate. It does not directly determine the subtraction magnitude.

Focal-cell intrinsic reference

Membrane markers

For membrane markers, the focal reference is the valid focal boundary that is not assigned to a plausible marker-positive source.

This makes broad membrane expression protective.

Example:

CD8 at source-facing interface: high
CD8 on the rest of focal membrane: also high
=> little interface excess
=> little or no correction

Conversely:

CD8 at source-facing interface: high
CD8 on the rest of focal membrane: near background
=> localized interface excess
=> contamination can be supported

Intracellular markers

For intracellular markers such as CD68, the self-reference comes from the eroded internal cell region.

Distributed internal signal therefore protects the cell, while a thin source-facing rim with little internal signal can be corrected.

Reference sufficiency

The package explicitly checks whether enough independent focal-cell reference remains.

Important settings:

"minimum_reference_valid_pixels": 4,
"minimum_reference_valid_fraction": 0.5,
"minimum_unconfounded_reference_fraction": 0.15,
"good_reference_fraction": 0.5

For membrane markers, minimum_unconfounded_reference_fraction is the important minimum: the focal cell needs enough valid membrane that is not occupied by plausible marker-positive source interfaces.

For intracellular markers, minimum_reference_valid_fraction governs the internal reference.

Reference quality is reported as values such as:

good

limited

insufficient

not_needed_no_plausible_source

marker_not_selected_for_correction

nuclear_reference_unavailable

Interface noise and evidence thresholds

The package estimates a robust marker-specific image noise scale and expresses interface evidence in units of that scale.

Default settings:

"interface_noise_threshold_floor_fraction": 0.05,
"interface_min_excess_noise_sd": 1.0,
"interface_strong_min_excess_noise_sd": 0.5,
"interface_high_specificity_min_excess_noise_sd": 2.0,
"interface_source_directionality_noise_sd": 1.0,
"interface_high_specificity_source_over_focal_noise_sd": 1.0

The standard interface model requires:

a plausible source;

a sufficient focal reference;

enough valid focal interface pixels;

focal interface excess above the standard noise threshold;

source signal above the focal reference by the configured directionality threshold.

The strong sensitivity method uses a more permissive interface-evidence threshold. The high_specificity method uses stricter interface evidence and also requires the source side to exceed the focal interface.

Physically supported contamination amount

For a supported directed pair, the core contamination amount is:

focal interface excess
x
valid focal interface pixel count / valid focal whole-cell pixel count

Therefore the amount removed comes from signal actually observed inside the focal segmentation.

A very bright neighboring cell cannot create an arbitrarily large subtraction if the focal interface contains only a small amount of excess signal.

Marker-attribution ambiguity

The CLI now has an explicit marker-level ambiguity output:

intrinsic_vs_neighbor_signal_ambiguous

It means:

The focal cell contains marker signal and plausible marker-positive neighbors occupy enough of the relevant reference region that the image cannot reliably distinguish intrinsic focal expression from neighbor-derived signal.

It does not mean that the cell is definitely double positive.

Important supporting fields include:

free_reference_fraction

source_contact_fraction

source_supported_signal_fraction

intrinsic_signal_support

n_plausible_source_neighbors

n_supported_interfaces

n_high_specificity_interfaces

reference_quality

ambiguity_reason

Default ambiguity settings:

"ambiguity_source_contact_fraction": 0.6,
"ambiguity_min_marker_positive_fraction": 0.05

Example: a CD8 T cell surrounded by CD4-positive cells can have unambiguous CD8 because its neighbors are not plausible CD8 sources, while CD4 on the same focal cell may be flagged ambiguous if CD4-positive neighbors occupy most of the usable membrane reference.

Ambiguous markers are preserved by the automatic recommendation.

Annotation policy

Correction is always annotation-free.

Supported modes are:

Mode

Behavior

disabled

Cell-type labels do not participate in correction

reporting_only

Labels may be carried into outputs for reporting

validation_only

Labels may be retained for downstream validation

weak_prior is no longer supported and will raise an error.

Even when celltype_col is configured, it does not alter pairwise interface correction or scenario recommendation.

Correction scenarios

All configured methods are calculated and saved. They are sensitivity views of the same pairwise-interface evidence, not seven independent biological models.

Scenario

Current meaning

none

No subtraction

conservative

Partial standard interface correction

medium

Full standard interface-supported correction

strong

Full correction using the more permissive supported-interface set

dominant_neighbor

Full standard correction from the largest supported source only

top_neighbors

Full standard correction from the top configured supported sources

high_specificity

Full correction using only high-specificity interfaces

Recommended scenario scaling:

"scenario_shrinkage": {
  "none": 0.0,
  "conservative": 0.5,
  "medium": 1.0,
  "strong": 1.0,
  "dominant_neighbor": 1.0,
  "top_neighbors": 1.0,
  "high_specificity": 1.0
}

The interface calculation is already physically bounded. Scenario maximum fractions therefore act only as emergency guardrails:

"scenario_max_fraction_removed": {
  "none": 0.0,
  "conservative": 1.0,
  "medium": 1.0,
  "strong": 1.0,
  "dominant_neighbor": 1.0,
  "top_neighbors": 1.0,
  "high_specificity": 1.0
}

Do not retain the old 0.25/0.50/0.80-style caps unless you deliberately want those arbitrary caps to override the new physically supported correction.

Automatic recommendation

The automatic recommendation is intentionally limited to:

none
conservative
medium

The other methods remain fully saved for manual selection and sensitivity analysis.

The decision rule is preservation-first:

marker not selected for correction
        -> none

intrinsic-vs-neighbor attribution ambiguous
        -> none + ambiguity flag

no supported contaminating interface
        -> none

supported contamination + substantial intrinsic focal signal
or limited but usable reference
        -> conservative

clear supported interface contamination + adequate reference
+ little intrinsic support
        -> medium

The marker-level threshold controlling whether intrinsic focal signal pushes the recommendation toward conservative is:

"recommendation_intrinsic_support_threshold": 0.25

strong, dominant_neighbor, top_neighbors, and high_specificity cannot automatically win recommendation through heuristic score bonuses.

All correction methods are retained

The authoritative long scenario table contains every configured scenario for every cell-protein pair:

correction_scenarios_<ROI>.parquet

This is intentional. A downstream QMD can compare:

Raw / none
Conservative
Medium
Strong
Dominant neighbor
Top neighbors
High specificity
Recommended

The recommended result does not replace or delete any alternative correction method.

Non-correction markers also remain present. Their corrected values are equal to their original values and they are labeled as not selected for correction.

Pair-level audit output

neighbor_contributions_<ROI>.parquet now represents directed pairwise interface evidence rather than whole-cell neighbor heuristics.

Depending on save_neighbor_contributions, it can contain fields such as:

focal label;

neighbor label;

protein;

focal interface intensity;

source interface intensity;

focal reference intensity;

valid interface pixels;

interface excess;

source plausibility;

standard/strong/high-specificity support flags;

physically supported contamination;

pair evidence strength;

source rank.

Configure storage with:

"save_neighbor_contributions": "top",
"max_saved_neighbors_per_cell_protein": 5

Use all only when full pair-level auditing is worth the larger output.

Geometry and density diagnostics

The workflow still calculates:

segmentation area;

number of touching neighbors;

total and maximum shared-boundary edge counts;

shared-boundary fraction proxy;

small-cell rank;

neighbor-density rank;

shared-boundary-density rank;

dense_small_cell_score.

These remain useful for QC. They no longer directly increase or decrease the correction amount.

Legacy keys such as dense_protection_strength are retained by the CLI only for compatibility with older configuration files and should not be treated as active correction controls.

Legacy correction keys

The current CLI still accepts several old keys so existing JSON files do not immediately fail, but they no longer control the pairwise-interface correction model:

annotation_prior_strength
minimum_neighbor_focal_contrast
strong_neighbor_focal_contrast
high_specificity_minimum_evidence
minimum_source_attribution_confidence
recommendation_minimum_margin
recommendation_minimum_confidence
allow_weighted_recommendation
dense_protection_strength
overcorrection_fraction_warning
retain_signed_corrected_values

A new project configuration should omit them unless compatibility with an older external wrapper requires them.

Signed and nonnegative corrected values are both written by the current implementation regardless of the legacy retain_signed_corrected_values switch.

Checkpoints and algorithm versioning

The correction algorithm version is included in checkpoint signatures beginning at stage 09.

Therefore, replacing an old installed CLI with the current pairwise-interface CLI automatically invalidates:

09_neighbor_exposure
10_correction_scenarios
11_recommendations
12_roi_h5ad
13_qc_plots
14_summary

while allowing valid stages 01-08 to be reused.

You normally do not need to set force_recompute_stages merely because the correction algorithm was updated.

Force a stage manually only when needed:

--force-stage 09_neighbor_exposure

For changes to registration, cropping, channel selection, or authoritative label mapping, force the appropriate earlier stage instead.

Main outputs

For ROI <ROI>, outputs are written under:

<outdir>/<ROI>/

Important outputs include:

Output

Description

spillover_features_<ROI>.parquet

Original per-cell protein/spatial features joined to metadata

cell_contact_pairs_<ROI>.parquet

Direct contact graph

geometry_density_features_<ROI>.parquet

Geometry and density diagnostics

neighbor_exposure_<ROI>.parquet

Cell-marker pairwise-interface exposure summary and ambiguity fields

neighbor_contributions_<ROI>.parquet

Optional directed source-interface audit records

correction_scenarios_<ROI>.parquet

Every saved correction scenario

suggested_corrections_<ROI>.parquet

Preservation-first recommendation and ambiguity metadata

correction_features_wide_<ROI>.parquet

Wide per-cell correction fields for AnnData

channel_thresholds_<ROI>.csv

Image-level marker thresholds

cell_id_to_raster_label_<ROI>.parquet

Table-label validation diagnostics

raster_mapping_summary_<ROI>.json

Mapping validation summary

roi_with_spillover_features_<ROI>.h5ad

ROI AnnData with original and correction outputs in .obs

spillover_summary_<ROI>.json

Final run summary

RUN_COMPLETE.json

Completion marker

qc_plots/

Image and feature QC

correction_qc/

Correction/recommendation QC

checkpoints/

Stage markers and cached artifacts

If Parquet support is unavailable, tables fall back to .csv.gz.

Wide AnnData correction fields

Scenario-specific columns are retained for all methods, for example:

protein_CD8A_none_corrected_nonnegative
protein_CD8A_conservative_corrected_nonnegative
protein_CD8A_medium_corrected_nonnegative
protein_CD8A_strong_corrected_nonnegative
protein_CD8A_dominant_neighbor_corrected_nonnegative
protein_CD8A_top_neighbors_corrected_nonnegative
protein_CD8A_high_specificity_corrected_nonnegative

The wide output also includes marker-level recommendation/ambiguity information such as:

protein_CD8A_suggested_corrected_nonnegative
protein_CD8A_suggested_fraction_removed
protein_CD8A_recommendation_confidence
protein_CD8A_n_plausible_source_neighbors
protein_CD8A_n_supported_interfaces
protein_CD8A_free_reference_fraction
protein_CD8A_source_contact_fraction
protein_CD8A_source_supported_signal_fraction
protein_CD8A_intrinsic_signal_support
protein_CD8A_intrinsic_vs_neighbor_signal_ambiguous

Do not discard original protein columns after correction.

QC interpretation

A well-behaved result should generally show:

no correction when there is no plausible marker-positive contacting source;

little correction when a membrane marker is broadly present around the focal cell;

stronger correction when excess signal is restricted to a source-facing interface;

correction bounded by the amount of focal-cell signal physically attributable to the interface;

preservation of cells lacking enough independent reference to resolve intrinsic versus neighbor-derived signal;

no systematic loss of true lineage markers merely because same-lineage cells cluster together;

all correction alternatives retained for comparison.

Important warnings include:

many cells with reference_quality = insufficient;

unexpectedly high intrinsic_vs_neighbor_signal_ambiguous rates;

large fractions removed from cells with strong intrinsic-signal support;

large corrections supported by very few valid interface pixels;

lineage-marker populations collapsing only under one aggressive sensitivity method;

systematic differences tied to ROI registration or QC-mask failure.

Troubleshooting

True lineage-positive cells are still overcorrected

Inspect:

free_reference_fraction;

intrinsic_signal_support;

focal interface means;

focal reference means;

n_supported_interfaces;

pair-level supported_contamination;

reference_quality;

the marker localization class.

Do not tune dense_protection_strength; it no longer controls correction.

Obvious contamination is undercorrected

Inspect:

source-interface positivity;

marker threshold;

robust noise scale;

interface excess;

source-over-reference directionality;

minimum_interface_valid_pixels;

interface_min_excess_noise_sd.

Too many ambiguity flags

Inspect whether cells genuinely lack unconfounded membrane/internal reference.

Potential settings are:

minimum_reference_valid_pixels
minimum_unconfounded_reference_fraction
minimum_reference_valid_fraction
ambiguity_source_contact_fraction
ambiguity_min_marker_positive_fraction

Do not lower these solely to eliminate flags. The ambiguity output exists specifically to represent cases the image cannot identify reliably.

Strong looks too aggressive

strong is a sensitivity output and is not eligible for automatic recommendation. Compare its pair-level support against medium rather than changing the recommendation system.

Neighbor-contribution tables are too large

Use:

"save_neighbor_contributions": "top",
"max_saved_neighbors_per_cell_protein": 5

or set save_neighbor_contributions to none.

Corrected signed values are negative

The signed output preserves the signed XOA measurement minus the estimated contamination. Use the separate nonnegative corrected field for downstream methods requiring nonnegative values.

Memory pressure

The new stage 09 performs pixel-level interface measurements, so it is more expensive than the old whole-cell exposure calculation. The implementation reuses memory-mapped image checkpoints and one ROI-level interface assignment rather than creating one image-sized mask per pair.

Practical controls remain:

analyze one ROI per job;

retain memory mapping;

use save_neighbor_contributions: top or none;

reduce unnecessary analysis_channels only if image-processing time becomes limiting;

increase Slurm memory/time only if real ROI profiling shows the current request is insufficient.

Downstream phenotype ambiguity review

The spillover CLI deliberately does not decide whether an ambiguous marker means a true double-positive or mixed-lineage cell.

After the CLI finishes, review_spillover_phenotype_ambiguity.py can be run as a separate post-processing step. It can compare phenotype assignments across all saved correction methods and prioritize cells whose correction uncertainty changes their biological interpretation.

This avoids manually reviewing every same-lineage cell in a crowded region. For example, a germinal-center B cell whose CD20 is neighbor-confounded but remains B-like under every correction method can be excluded from review, whereas a T cell switching between CD8-only and CD4/CD8-double-positive can be prioritized.

Recommended project workflow

Verify registration and authoritative table-to-raster mapping.

Run one representative ROI.

Confirm XOA preprocessing and QC-mask behavior.

Inspect original cell-level protein features.

Inspect pairwise interface evidence for CD3E, CD4, CD8A, CD20, CD68, and other key lineage markers.

Confirm same-marker neighboring cells are preserved.

Confirm clear cross-lineage interface contamination is corrected.

Inspect ambiguity flags in densely packed immune regions.

Compare all seven correction methods in the QMD; do not hide sensitivity outputs behind recommended.

Confirm the recommended result uses only none, conservative, or medium.

Expand to all ROIs.

Run the downstream phenotype-ambiguity reviewer on completed outputs.

Perform targeted manual review only for biologically consequential ambiguous cells.

Full option reference

protein-spillover run --help

Generate the exact current configuration template with:

protein-spillover template-config \
  --output full_protein_spillover_config.json