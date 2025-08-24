import subprocess
import time

sra_numbers = [
    "SRR7179504", "SRR7179505", "SRR7179506", "SRR7179507",
    "SRR7179508", "SRR7179509", "SRR7179510", "SRR7179511",
    "SRR7179520", "SRR7179521", "SRR7179522", "SRR7179523",
    "SRR7179524", "SRR7179525", "SRR7179526", "SRR7179527", "SRR7179536", "SRR7179537", "SRR7179540","SRR7179541"
    ]

for sra_id in sra_numbers:
    print("\n=== Downloading:", sra_id, "===")
    prefetch_cmd = f"prefetch {sra_id}"
    print("Command:", prefetch_cmd)

    start_time = time.time()
    subprocess.call(prefetch_cmd, shell=True)
    end_time = time.time()

    elapsed_min = (end_time - start_time) / 60
    print(f"⏱ Download time for {sra_id}: {elapsed_min:.2f} minutes")
for sra_id in sra_numbers:
    sra_path = f"./{sra_id}/{sra_id}.sra"
    print("\n=== Generating FASTQ for:", sra_id, "===")

    # Using fasterq-dump (faster + multithreaded)
    fastq_dump_cmd = (
        f"fasterq-dump {sra_path} --split-files -O fastq --threads 8"
    )
    print("Command:", fastq_dump_cmd)

    start_time = time.time()
    subprocess.call(fastq_dump_cmd, shell=True)

    # Compress output FASTQs with pigz (parallel gzip)
    compress_cmd = f"pigz -p 8 fastq/{sra_id}*.fastq"
    subprocess.call(compress_cmd, shell=True)

    end_time = time.time()
    elapsed_min = (end_time - start_time) / 60
    print(f"⏱ FASTQ generation time for {sra_id}: {elapsed_min:.2f} minutes")