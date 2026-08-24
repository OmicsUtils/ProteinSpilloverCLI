# Spatial protein spillover CLI

`protein-spillover` is an end-to-end workflow for measuring, diagnosing, and correcting cell-segmentation-associated protein spillover in Xenium and other multiplex protein images stored in `SpatialData`.

The workflow does more than report one corrected intensity. For every cell and protein, it preserves the original measurement, extracts spatial evidence for possible spillover, evaluates multiple correction scenarios, estimates overcorrection risk, and recommends a corrected value with structured and human-readable reasons.

The default correction path is annotation-free. Cell-type labels are optional and are never required.

## What the workflow does

The CLI performs the complete pipeline in one program:

1. loads the SpatialData table, protein image, label raster, and optional QC mask;
2. validates the authoritative table-to-raster label mapping;
3. crops the selected ROI on the verified native image grid;
4. preprocesses image intensities correctly for generic, precorrected, or Xenium XOA inputs;
5. extracts whole-cell, interior, boundary, outer-ring, and neighboring-cell protein features;
6. builds a direct cell-contact graph;
7. calculates segmentation geometry, local density, and dense-small-cell protection scores;
8. estimates marker-specific neighbor exposure and likely spillover sources;
9. evaluates multiple correction scenarios;
10. recommends a scenario or marks the result unresolved;
11. saves confidence, overcorrection risk, reason codes, and readable explanations;
12. writes corrected features into an ROI AnnData object and produces QC reports.

The input Zarr remains read-only. Results are written only to the configured output directory.

## Core assumptions

The workflow is appropriate when all of the following are true:

1. The table has one row per segmented biological cell.
2. A positive integer table column directly identifies the same cell in the label raster.
3. The protein image and label raster share the same full-resolution native `y × x` grid.
4. Table spatial coordinates can be related to native image pixels using an axis-aligned transformation consisting of:
   - one pixel-size scale;
   - optional x and/or y reflection;
   - optional global x and y origins.
5. The original image values should remain available for auditing.
6. The workflow should not write back into the authoritative SpatialData Zarr.

The CLI does not infer a label crosswalk from centroid matching or polygon overlap. The table label column is treated as authoritative and checked against the cropped label raster.

## Files in this distribution

```text
protein_spillover_cli.py       Main end-to-end CLI and analysis workflow
protein_spillover_config.json  Complete configuration template
protein_spillover_README.md    This guide
requirements.txt               Dependency list
pyproject.toml                 Editable installation and console entry point
```

## Installation

Use Python 3.10 or newer in an isolated environment containing the SpatialData stack required by the project.

### Install as a command

From the distribution directory:

```bash
python -m pip install -e .
```

This installs:

```bash
protein-spillover
```

### Run without installing

```bash
python protein_spillover_cli.py --help
```

### Verify installation

```bash
protein-spillover --help
protein-spillover run --help
```

## First inspect a new project

Use `inspect` before building a production configuration.

```bash
protein-spillover inspect \
  --zarr /data/project/spatial_data.zarr
```

A more detailed inspection:

```bash
protein-spillover inspect \
  --zarr /data/project/spatial_data.zarr \
  --table cell_table \
  --image morphology_focus \
  --labels cell_labels \
  --roi-col roi
```

Use this output to identify:

- the table name;
- protein image and label element names;
- the stable cell-ID column;
- the authoritative integer raster-label column;
- the `obsm` key containing x/y coordinates;
- available protein channels;
- optional ROI, sample, batch, and annotation columns.

## Create a configuration file

Generate a complete template:

```bash
protein-spillover template-config \
  --output protein_spillover_config.json
```

Run with the configuration:

```bash
protein-spillover run \
  --config protein_spillover_config.json
```

CLI values override JSON values:

```bash
protein-spillover run \
  --config protein_spillover_config.json \
  --roi-value ROI_002 \
  --outdir /results/protein_spillover/ROI_002
```

## Minimal direct-argument run

```bash
protein-spillover run \
  --zarr /data/project/spatial_data.zarr \
  --table cell_table \
  --image morphology_focus \
  --labels cell_labels \
  --outdir /results/protein_spillover \
  --cell-id-col cell_id \
  --table-label-col cell_labels \
  --spatial-key spatial \
  --pixel-size-um 0.2125 \
  --orientation no_flip \
  --roi-col roi \
  --roi-value ROI_001 \
  --intensity-mode xenium_xoa \
  --channel CD3E \
  --channel CD8A
```

To analyze the full table rather than one ROI:

```bash
protein-spillover run \
  --config protein_spillover_config.json \
  --all-cells
```

When `roi_col` is configured but no ROI value is supplied, the workflow chooses the ROI whose cell count is closest to the median ROI cell count. This is intended for pilot testing, not for automatic all-ROI processing.

## Required configuration

| JSON key | CLI flag | Meaning |
|---|---|---|
| `sdata_zarr_path` | `--zarr` | Input SpatialData Zarr |
| `table_name` | `--table` | AnnData table stored in SpatialData |
| `protein_image_name` | `--image` | Multichannel protein image |
| `cell_labels_name` | `--labels` | Positive-integer cell-label raster |
| `outdir` | `--outdir` | Output directory |
| `cell_id_col` | `--cell-id-col` | Stable cell identifier in `table.obs` |
| `table_cell_label_col` | `--table-label-col` | Authoritative raster label in `table.obs` |
| `spatial_key` | `--spatial-key` | `table.obsm` key containing x/y coordinates |
| `native_pixel_size_um` | `--pixel-size-um` | Native pixel size in table coordinate units |

Table label values must be unique positive integers within the selected ROI.

## Image registration

The workflow uses an explicit axis-aligned native-pixel mapping rather than trusting stored transformations.

### Orientation

| Value | Meaning |
|---|---|
| `no_flip` | x and y increase with native image indices |
| `x_flip` | x is reflected |
| `y_flip` | y is reflected |
| `xy_flip` | both axes are reflected |

Arbitrary rotation, shear, and nonlinear warping are not supported. Such images should be registered before this workflow runs.

### Pixel size and origins

For `no_flip`:

```text
x_table = native_origin_x_um + x_native × native_pixel_size_um
y_table = native_origin_y_um + y_native × native_pixel_size_um
```

Relevant options:

```text
--pixel-size-um
--origin-x-um
--origin-y-um
--orientation
```

Before a production run, verify registration on several spatially separated cells.

## ROI behavior

### Named ROI

```json
{
  "roi_col": "roi",
  "pilot_roi": "ROI_001"
}
```

### Representative ROI

Set `roi_col` and leave `pilot_roi` as `null`.

### Full table

```json
{
  "roi_col": null,
  "pilot_roi": null
}
```

or pass `--all-cells`.

## Channel selection

Select individual channels:

```bash
--channel CD3E --channel CD8A --channel CD20
```

Use all channels except exclusions:

```json
"analysis_channels": null,
"exclude_channels": ["DAPI"]
```

Repeated CLI channel options replace the JSON list.

## Intensity preprocessing

The workflow supports three modes.

### `xenium_xoa`

Use this for Xenium Onboard Analysis protein images.

XOA protein images are already deconvolved, autofluorescence-background-subtracted, saturation-masked, and spectrally corrected. XOA normally stores these values with a positive offset of 100.

The workflow:

- subtracts the configured XOA storage offset;
- retains signed offset-adjusted values for intensity summaries;
- creates a clipped nonnegative signal representation for thresholds and ratios;
- does not apply a second Gaussian background subtraction.

```json
{
  "input_intensity_mode": "xenium_xoa",
  "xenium_xoa_intensity_offset": 100.0
}
```

### `precorrected`

Use supplied values directly as signed analysis intensities.

### `generic_gaussian`

Estimate and subtract a broad Gaussian background for generic multiplex images that have not already been background corrected.

## QC masks

A channel-matched QC-mask image may be supplied.

For Xenium XOA saturation masks, valid pixels are usually 0 and masked pixels are 255.

```json
{
  "protein_qc_mask_name": "morphology_focus_qc_masks",
  "qc_mask_valid_value": 0,
  "xenium_require_qc_mask": false
}
```

Invalid pixels are excluded from intensity summaries, thresholds, directionality metrics, and QC distributions.

## Cell masks and spatial compartments

For every segmented cell, the workflow derives:

- whole cell;
- eroded interior;
- internal boundary;
- outer ring;
- outer noncell ring;
- outer other-cell ring.

Configure with:

```bash
--inner-erosion 2
--outer-ring 2
```

These values are full-resolution native pixels.

## Exploratory positive-pixel thresholds

These thresholds are used for image-localization and directionality evidence. They are not biological protein gates.

Default quantile:

```bash
--threshold-quantile 0.90
```

Channel-specific quantiles:

```bash
--channel-threshold-quantile CD3E=0.85 \
--channel-threshold-quantile PD1=0.95
```

Manual image thresholds:

```bash
--manual-threshold CD3E=120
```

Manual thresholds take precedence.

## Annotation policy

Cell-type labels are optional.

The default is:

```json
"annotation_mode": "disabled"
```

Supported modes:

| Mode | Behavior |
|---|---|
| `disabled` | Labels do not affect correction or recommendation |
| `reporting_only` | Labels are copied to outputs only |
| `validation_only` | Labels are used only in QC summaries |
| `weak_prior` | Labels may weakly influence recommendation scoring |

The annotation-free result remains the primary technical correction path. Labels should never force a marker to zero or define a correction by themselves.

## Geometry and dense-small-cell protection

The workflow computes per-cell and local-neighborhood features including:

- segmentation area;
- perimeter and shared-boundary exposure;
- number of direct contacts;
- total and maximum shared-boundary edge counts;
- local neighbor density;
- small-cell rank;
- dense-neighborhood rank;
- dense-small-cell geometry score.

This score is deliberately geometric. It does not claim that a cell is a lymphocyte.

Dense-small-cell protection reduces correction aggressiveness when a cell is small, crowded, and surrounded by several similarly exposed neighbors. This is intended to limit density-dependent overcorrection in lymphocyte-rich regions.

Main settings:

```json
{
  "dense_small_cell_area_quantile": 0.35,
  "dense_neighbor_count_quantile": 0.70,
  "dense_shared_boundary_quantile": 0.70,
  "dense_protection_strength": 0.65
}
```

## Neighbor exposure and source attribution

For each cell-protein pair, the workflow evaluates directly contacting cells using:

- shared-boundary evidence;
- neighbor intensity;
- neighbor-to-focal contrast;
- dominant-neighbor contribution;
- top-neighbor contributions;
- boundary-versus-interior signal localization;
- anisotropy and angular boundary coverage;
- dense-small-cell protection.

The workflow distinguishes:

- one dominant bright source;
- several strong sources;
- many weak neighbors;
- dense homogeneous neighborhoods with no clear source.

This reusable evidence is calculated once and reused by all correction scenarios.

## Correction scenarios

The workflow supports a configurable scenario engine. At minimum, the required anchors are `none`, `conservative`, `medium`, and `strong`.

Default scenarios:

| Scenario | Purpose |
|---|---|
| `none` | Preserve the original value when spillover evidence is weak |
| `conservative` | Protect dense small-cell regions and limit removed signal |
| `medium` | Standard correction for moderate directional evidence |
| `strong` | Aggressive correction for compelling, high-confidence spillover |
| `dominant_neighbor` | Use only the strongest plausible source |
| `top_neighbors` | Use only the strongest configured number of sources |
| `high_specificity` | Correct only when evidence exceeds a stringent threshold |

Scenario behavior is controlled by:

```json
"scenario_shrinkage": {
  "none": 0.0,
  "conservative": 0.30,
  "medium": 0.60,
  "strong": 0.90,
  "dominant_neighbor": 0.70,
  "top_neighbors": 0.65,
  "high_specificity": 0.75
}
```

and by maximum removable fractions:

```json
"scenario_max_fraction_removed": {
  "none": 0.0,
  "conservative": 0.25,
  "medium": 0.50,
  "strong": 0.80,
  "dominant_neighbor": 0.55,
  "top_neighbors": 0.60,
  "high_specificity": 0.65
}
```

These are sensitivity scenarios, not claims that one universal correction strength is correct for every cell.

## Suggested corrected value

For each cell-protein pair, the workflow selects one applicable scenario or returns `unresolved`.

The recommendation considers:

- evidence that spillover exists;
- source-attribution confidence;
- segmentation confidence;
- boundary localization;
- neighbor-to-focal contrast;
- dense-small-cell protection;
- proposed fraction removed;
- disagreement between scenarios;
- recommendation margin between the best and second-best scenarios.

The output includes:

- suggested scenario;
- suggested signed corrected intensity;
- suggested nonnegative corrected intensity;
- estimated contamination removed;
- second-best scenario;
- selection margin;
- recommendation confidence;
- overcorrection risk;
- structured reason codes;
- readable reason text.

The suggested value remains auditable because all alternative scenarios are retained.

## Confidence fields

Confidence is intentionally decomposed.

| Field | Meaning |
|---|---|
| segmentation confidence | Reliability of the focal cell mask and derived geometry |
| localization confidence | Strength of boundary-versus-interior evidence |
| source-attribution confidence | Ability to assign spillover to one or more neighbors |
| bleeding confidence | Evidence that meaningful spillover is present |
| recommendation confidence | Confidence in selecting one correction scenario |
| overcorrection risk | Risk that subtraction removes genuine focal-cell signal |

These values do not represent biological protein-positivity confidence. Biological gating belongs in the downstream gating and phenotyping workflow.

## Neighbor-contribution output

Per-neighbor output may become large.

```json
{
  "save_neighbor_contributions": "top",
  "max_saved_neighbors_per_cell_protein": 5
}
```

Modes:

| Value | Behavior |
|---|---|
| `none` | Save aggregate exposure only |
| `top` | Save the strongest contributors per cell-protein pair |
| `all` | Save every eligible neighbor contribution |

Use `all` cautiously on large datasets.

## Checkpoints and restarting

The workflow contains fourteen ordered stages:

```text
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
```

Resume is enabled by default. Each checkpoint includes a configuration signature and upstream signature.

Changing correction parameters invalidates only the correction and downstream stages. It does not force image cropping or pixel-feature extraction to rerun.

Disable resume:

```bash
--no-resume
```

Force one stage and every downstream stage:

```bash
--force-stage 10_correction_scenarios
```

Large array checkpoints:

```bash
--save-array-checkpoints
--memory-map-arrays
--cleanup-arrays
```

## Main outputs

For ROI `<ROI>`, outputs are written beneath:

```text
<outdir>/<ROI>/
```

Important outputs include:

| Output | Description |
|---|---|
| `spillover_features_<ROI>.parquet` | Original per-cell protein and spatial features joined to metadata |
| `cell_contact_pairs_<ROI>.parquet` | Direct contact graph with shared-boundary evidence |
| `geometry_density_<ROI>.parquet` | Cell geometry, contact density, and dense-small-cell scores |
| `neighbor_exposure_<ROI>.parquet` | Cell-protein spillover evidence and source attribution |
| `neighbor_contributions_<ROI>.parquet` | Optional retained source-neighbor contributions |
| `correction_scenarios_<ROI>.parquet` | All scenario-level corrected values and contamination estimates |
| `suggested_corrections_<ROI>.parquet` | Recommended scenarios, confidence, risks, and reasons |
| `correction_features_wide_<ROI>.parquet` | Per-cell wide correction columns for AnnData integration |
| `channel_thresholds_<ROI>.csv` | Exploratory image thresholds by marker |
| `cell_id_to_raster_label_<ROI>.parquet` | Table-label validation and centroid diagnostics |
| `raster_mapping_summary_<ROI>.json` | Mapping coverage and validation summary |
| `roi_with_spillover_features_<ROI>.h5ad` | ROI AnnData with original features and correction outputs in `.obs` |
| `spillover_summary_<ROI>.json` | Final run summary |
| `RUN_COMPLETE.json` | Successful completion marker |
| `qc_plots/` | Intensity and feature QC |
| `correction_qc/` | Correction magnitude, density, and recommendation QC |
| `checkpoints/` | Stage markers and cached artifacts |

If Parquet support is unavailable, tables fall back to compressed CSV files ending in `.csv.gz`.

## Feature naming

Original image-derived features use names such as:

```text
protein_<channel>_raw_whole_mean
protein_<channel>_whole_mean
protein_<channel>_inner_mean
protein_<channel>_boundary_mean
protein_<channel>_outer_othercell_mean
protein_<channel>_boundary_to_inner_log2
protein_<channel>_boundary_angular_coverage
protein_<channel>_boundary_anisotropy
```

Correction features written to AnnData include scenario-specific and suggested outputs using safe marker and scenario names.

Do not discard the original columns after correction.

## QC interpretation

A well-behaved correction should generally show:

- limited correction when no brighter neighbor is present;
- stronger correction when one bright neighbor dominates and boundary evidence is directional;
- reduced aggressiveness in dense small-cell regions;
- no severe downward trend in corrected lineage-marker intensity with local density;
- preservation of isolated high-intensity cells;
- a manageable unresolved fraction;
- few extreme negative signed corrected values.

Warnings include:

- corrected intensity decreasing sharply as local density increases;
- strong corrections driven by many weak neighbors;
- large fractions removed without a dominant source;
- widespread strong correction in dense immune aggregates;
- scenario disagreement across a large fraction of cells;
- low source-attribution confidence with high recommended correction.

## Common troubleshooting

### Table labels are absent from the crop

Check pixel size, orientation, origins, ROI coordinates, table-label column, crop margin, and whether the selected image and labels belong to the same acquisition.

Do not lower `minimum_raster_mapping_fraction` merely to bypass registration problems.

### Dense lymphocytes appear overcorrected

Review:

- dense-small-cell score;
- neighbor count;
- shared-boundary fraction;
- dominant-source fraction;
- correction scenario selected;
- fraction removed;
- correction versus density QC plots.

Consider increasing `dense_protection_strength`, lowering conservative shrinkage, or reducing the conservative maximum removable fraction. Re-run from `08_geometry_density` or `10_correction_scenarios`, depending on which settings changed.

### Strong correction is selected too often

Increase `strong_neighbor_focal_contrast`, increase `recommendation_minimum_confidence`, or increase `recommendation_minimum_margin`.

### Too many unresolved results

Inspect whether scenario scores are genuinely close. Lowering the minimum margin may reduce unresolved calls but can create false precision.

### Neighbor-contribution tables are too large

Use `save_neighbor_contributions: top` or `none` and reduce `max_saved_neighbors_per_cell_protein`.

### Corrected values are negative

Signed corrected values are retained intentionally. Use the separate nonnegative corrected output when a downstream method requires nonnegative values. Do not silently overwrite signed values.

### Memory pressure

Analyze one ROI at a time, reduce channels, retain memory mapping, use top-neighbor output, or increase QC downsampling. Scenario evaluation itself should be relatively inexpensive because it reuses cached exposure evidence.

## Recommended workflow for a new project

1. Run `inspect` and record exact element, table, column, and channel names.
2. Verify the authoritative table label against the raster.
3. Confirm image-label native-grid equality.
4. Verify pixel size, orientation, and origins on several cells.
5. Start with one representative ROI and a small channel set.
6. Review raster-mapping and intensity QC.
7. Review whole-cell, interior, boundary, and outer-other-cell features.
8. Review geometry-density and neighbor-exposure outputs.
9. Compare all correction scenarios before trusting the suggested value.
10. Specifically examine dense lymphocyte-rich regions for overcorrection.
11. Only then expand to all channels and ROIs.
12. Pass original, alternative, and suggested corrected values to the downstream gating and phenotyping CLI.

## Full option reference

```bash
protein-spillover run --help
```

Generate the exact configuration schema:

```bash
protein-spillover template-config \
  --output full_protein_spillover_config.json
```
