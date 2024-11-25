"""
Simple script to download
"""

import os
import argparse
import multiprocessing

def download_from_directory(task_id, input_directory):
    """
    Function to be executed in parallel that downloads all the files for a specific run
    """

    if "hy_" not in input_directory:
        return

    print(f"Processing task {task_id}")
    train_id = input_directory.split(sep="/")[-1]
    os.system(f"alien_find alien://{input_directory} "
              f"AOD/*/AnalysisResults.root > outputs_{train_id}.txt")

    check_unmerged = False
    with open(f"outputs_{train_id}.txt") as file:
        lines = [line.rstrip() for line in file]

        if len(lines) == 0:
            check_unmerged = True

    if check_unmerged:
        os.system(f"alien_find alien://{input_directory} "
                  f"*/AnalysisResults.root > outputs_{train_id}.txt")
        with open(f"outputs_{train_id}.txt") as file:
            lines = [line.rstrip() for line in file]

    if not os.path.isdir(train_id):
        os.mkdir(train_id)
    for ifile, line in enumerate(lines):
        if "AnalysisResults" not in line:
            continue
        os.system(f"alien_cp {line} file:{train_id}/AnalysisResults_{ifile}.root")
    print(f"Processing task {task_id} - DONE")


def download_and_merge(input_file, num_workers, output_file):
    """
    Main function to download and merge output files from hyperloop
    """

    output_directories = []
    with open(input_file, 'r') as file:
        for line in file:
            # Split line by comma or any delimiter
            parts = line.strip().split(',')
            for part in parts:
                output_directories.append(part)

    num_workers = min(num_workers, os.cpu_count())  # Get the number of available CPU cores)

    # tasks = [(idir, path) for (idir, path) in enumerate(output_directories)]
    with multiprocessing.Pool(processes=num_workers) as pool:
        pool.starmap(download_from_directory, enumerate(output_directories))

    if len(output_file.split("/")) != 1:
        prev_path = ""
        for folder in output_file.split("/")[:-1]:
            os.mkdir(f"{prev_path}{folder}")
            prev_path += f"{folder}/"

    os.system(f"hadd -f {output_file} */AnalysisResults*.root")
    os.system(f"rm -r hy_*")
    os.system(f"rm outputs_hy_*")

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description="Arguments")
    PARSER.add_argument("--input_file", "-i", metavar="text",
                        help="text input file with directories", required=True)
    PARSER.add_argument("--jobs", "-j", type=int, default=20,
                        help="number of workers", required=False)
    PARSER.add_argument("--output_file", "-o", metavar="text",
                        default="AnalysisResults_LHC23_Dstar_train293534.root",
                        help="output file", required=False)
    ARGS = PARSER.parse_args()

    download_and_merge(ARGS.input_file, ARGS.jobs, ARGS.output_file)
