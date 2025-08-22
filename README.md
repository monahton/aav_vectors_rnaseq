

#  README

## 1. Data sources

This task is based on publicly available sequencing data from a study of **\[insert study topic here]**. The dataset includes multiple samples under different conditions (e.g., treated vs. control) and was originally sequenced using **\[insert platform, e.g., Illumina 2×150]**.
The subsampled and cleaned FASTQs are stored in `data/` and are used as the inputs for the workflow.

---

## 2. How to download

INSTRUCTIONS TO ACCESS THE DATA
### Example using SRA Toolkit

```bash
CODE TO DOWNLOAD
```


---

## 3. Pre-processing / subsampling

INCLUDE THE METHOD YOU USED TO SUBSAMPLE, MINATURIZE, OR TRIM DOWN

1. **STEP 1** ...

Example:

```bash
CODE TO SUBSAMPLE
```


---

## 4. How the workflow works
DESCRIBE THE WORKFLOW HERE - NOTE THE BELOW ARE JUST EXAMPLES, REPLACE WITH YOUR OWN - YOURS CAN TAKE A VERY DIFFERENT FORMAT
The workflow files is stored in `workflow/`.

---

### Step 1 – Quality Control (example)

**Purpose:** Remove low-quality reads and adapter sequences
**Tools:** `fastp`, `cutadapt`, `trimmomatic`
**Inputs:** Subsampled FASTQ files (from `data/fastq_subsampled/`)
**Outputs:** Cleaned FASTQs, QC reports (`.html`, `.json`, or `.txt`)
**Command:**

```bash
fastp --in1 sample.fastq.gz --out1 cleaned.fastq.gz ...
```

---

### Step 2 ...

**Purpose:** ...
**Tools:** ...
**Inputs:** ...
**Outputs:** ...
**Command:**


---

### Step X – Analysis (e.g., DESeq2, variant calling, etc.)

**Purpose:** ...
**Tools:** ...
**Inputs:** ...
**Outputs:** ...
**Command:**

