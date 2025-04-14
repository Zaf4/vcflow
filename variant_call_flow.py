import zipfile
from pathlib import Path
from subprocess import run

import requests
from tqdm import tqdm

"""
TODO:
1. unzip the zip file
2. get the reference genome hg38 from UCSC
3. Index the reference genome (bwa index)
3. align the reads to the reference genome (bwa mem)
4. post-process the alignment
5. call variants (GATK)
"""


def red(text: str) -> str:
    """Return red text."""
    return f"\033[31m{text}\033[0m"


def cyan(text: str) -> str:
    """Return cyan text."""
    return f"\033[36m{text}\033[0m"


def unzip(file_name: str) -> None:
    """Unzip a zip file."""
    target_dir = Path(file_name.replace(".zip", ""))  # remove the .zip extension
    if target_dir.exists():
        print(f"{target_dir} already exists, skipping unzipping...")
    else:
        with zipfile.ZipFile(file_name, "r") as zip_file:
            zip_file.extractall()

    return


def get_reference_genome(url: str, gunzip: bool = True) -> None:
    """Get the reference genome hg38 from UCSC."""
    ref_dir = Path("ref")
    ref_dir.mkdir(exist_ok=True)
    response = requests.get(url, stream=True)
    output_path = ref_dir / url.split("/")[-1]  # ref/last part of the url
    # Check if the reference genome file already exists
    if output_path.exists():
        print(f"{output_path} already exists, skipping downloading...")
    else:
        total_size = int(response.headers.get("content-length", 0))
        # Download the reference genome file
        with (
            Path.open(output_path, "wb") as f,
            tqdm(  # progress bar
                desc="Downloading",
                total=total_size,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
            ) as bar,
        ):
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    bar.update(len(chunk))
    if gunzip:
        run(["gunzip", output_path, "--keep"])

    return


def index_reference_genome(ref: str) -> None:
    """
    Index the reference genome.

    output will be:

    *..amb	Ambiguity file: stores ambiguous bases (like N)
    *..ann	Annotation file: contains names and lengths of sequences
    *..bwt	BWT-transformed sequence used for fast alignment
    *..pac	Packed FASTA: 4 bases per byte; helps with space efficiency
    *..sa	Suffix array: enables fast substring matching
    """
    ref_path = Path(ref)
    suffix_array = Path(ref + ".sa")
    # if the index files already exist, skip indexing
    if suffix_array.exists():  # the last file is .sa which is the suffix array
        print(f"{suffix_array} already exists, skipping indexing...")
    else:
        print(cyan(f"indexing {ref} with bwa"))
        run(["bwa", "index", ref_path])

    return


def faidx(ref: str) -> None:
    """Index the reference genome with samtools."""
    ref_path = Path(ref)
    fai_path = Path(ref + ".fai")
    # if the index files already exist, skip indexing
    if fai_path.exists():
        print(f"{fai_path} already exists, skipping indexing...")
    else:
        run(["samtools", "faidx", ref_path])


def align_reads_paired_end(
    read1: str, read2: str, ref: str, output: str, n_threads: int = 8
) -> None:
    """Align reads to the reference genome."""
    ref_path = Path(ref)
    output_path = Path(output)
    output_path.parent.mkdir(exist_ok=True, parents=True)
    # if the alignment files already exist, skip alignment
    if output_path.exists():
        print(f"{output_path} already exists, skipping alignment...")
    else:
        with Path.open(output_path, "w") as output_file:
            run(
                ["bwa", "mem", "-t", str(n_threads), ref_path, read1, read2],
                stdout=output_file,
            )

    return


def postprocess_alignment(sam_file: str) -> None:
    """Post process the alignment."""
    sam_file_path = Path(sam_file)
    bam_file_path = sam_file_path.with_suffix(".bam") # .bam
    read_grouped_bam_file_path = bam_file_path.with_suffix(".rg.bam") # .rg.bam
    sorted_bam_file_path = read_grouped_bam_file_path.with_suffix(".sorted.bam") # .rg.sorted.bam
    stats_file_path = sorted_bam_file_path.with_suffix(".bam.stats") # .rg.sorted.bam.stats
    index_file_path = sorted_bam_file_path.with_suffix(".bam.bai") # .rg.sorted.bam.bai

    # covert sam to bam
    if bam_file_path.exists():
        print(f"{bam_file_path} already exists, skipping conversion...")
    else:
        print(cyan(f"converting {sam_file_path} to {bam_file_path}"))
        with Path.open(bam_file_path, "w") as bam_file:
            run(["samtools", "view", "-bhS", sam_file_path], stdout=bam_file)
    # add read groups
    if not read_grouped_bam_file_path.exists():
        print(cyan(f"adding read groups to {bam_file_path}"))
        run(
            [
                "samtools",
                "addreplacerg",
                bam_file_path,
                "-r",
                "@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA",  # customize these RG tags
                "-o",
                read_grouped_bam_file_path,
            ]
        )
    # sort the bam file
    if sorted_bam_file_path.exists():
        print(f"{sorted_bam_file_path} already exists, skipping sorting...")
    else:
        print(cyan(f"sorting {bam_file_path} to {sorted_bam_file_path}"))
        run(["samtools", "sort", read_grouped_bam_file_path, "-o", sorted_bam_file_path])
    # generate stats
    if stats_file_path.exists():
        print(f"{stats_file_path} already exists, skipping stats...")
    else:
        print(cyan(f"generating stats for {sorted_bam_file_path}"))
        with Path.open(stats_file_path, "w") as stats_file:
            run(["samtools", "stats", sorted_bam_file_path], stdout=stats_file)
    # index the bam file
    if index_file_path.exists():
        print(f"{index_file_path} already exists, skipping indexing...")
    else:
        print(cyan(f"indexing {sorted_bam_file_path}"))
        run(["samtools", "index", sorted_bam_file_path])

    return


def call_variants(bam_file: str, ref: str, regions_bed: str, output_vcf: str) -> None:
    """Call variants."""
    vcf_file_path = Path(output_vcf)
    ref_path = Path(ref)

    dict_file_path = ref_path.with_suffix(".dict")
    # if the vcf file already exists, skip calling variants
    if not dict_file_path.exists():
        print(cyan(f"creating dict file for {ref}"))
        run(["gatk", "CreateSequenceDictionary", "-R", ref, "-O", dict_file_path])
    if vcf_file_path.exists():
        print(f"{vcf_file_path} already exists, skipping calling variants...")
    else:
        print(cyan(f"calling variants in {bam_file} with {regions_bed}"))
        # finally call variants
        run(
            [
                "gatk",
                "HaplotypeCaller",
                "-R",
                ref_path,
                "-I",
                bam_file,
                "-L",
                regions_bed,
                "-O",
                output_vcf,
            ]
        )

    return


def main() -> None:
    """Main function to combine the flow."""
    unzip(file_name="Task_Files.zip")
    get_reference_genome(
        url="https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz", gunzip=True
    )
    index_reference_genome(ref="ref/hg38.fa")
    align_reads_paired_end(
        read1="Task_Files/POOL-37_S182_L004_R1_001.fastq.gz",
        read2="Task_Files/POOL-37_S182_L004_R2_001.fastq.gz",
        ref="ref/hg38.fa",
        output="output/POOL-37.sam",
    )
    postprocess_alignment(sam_file="output/POOL-37.sam")
    call_variants(
        bam_file="output/POOL-37.sorted.bam",
        ref="ref/hg38.fa",
        regions_bed="Task_Files/POOL-37.bed",
        output_vcf="output/POOL-37.vcf",
    )
    return


if __name__ == "__main__":
    main()
