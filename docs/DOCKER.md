# Docker documentation

If you have Docker installed, pull the Docker image from Docker Hub.

```bash
docker pull sfilges/umierrorcorrect2
```

Download a reference genome fasta file and mount the reference directory and data directory (including fastq files and BED files) to the docker container. The container is configured with `umierrorcorrect` as the entrypoint.

To view the help message:

```bash
docker run -it sfilges/umierrorcorrect2 --help
```

## Running the pipeline

Since the umierrorcorrect pipeline uses `bwa` for mapping reads, a bwa-indexed reference genome is needed.

To run the full pipeline, use the `batch` command. Note that all file paths must be relative to the mapped volumes inside the container (e.g., `/data/...` or `/references/...`).

Example command:

```bash
docker run -v /path/to/references/:/references/ -v /path/to/data/:/data/ -it sfilges/umierrorcorrect2 \
    batch \
    -r1 /data/read1.fastq.gz \
    -r2 /data/read2.fastq.gz \
    -ul 12 \
    -sl 0 \
    -r /references/reference.fa \
    -o /data/output_directory
```

You can also run individual steps, for example `preprocess` or `consensus` instead of `batch`.
