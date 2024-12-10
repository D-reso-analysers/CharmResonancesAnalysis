import ROOT
import os
import argparse

# Function to copy a directory from one ROOT file to another
def remove_directory(source_file, target_file, dir_name):
    # Open the source file
    src = ROOT.TFile.Open(source_file, "UPDATE")
    if not src or src.IsZombie():
        raise RuntimeError(f"Cannot open source file: {source_file}")
    
    # Check if the directory exists
    src_dir = src.Get(dir_name)
    if not src_dir:
        raise RuntimeError(f"Directory '{dir_name}' not found in {source_file}")
    
    # Open the target file
    tgt = ROOT.TFile(target_file, "recreate")
    

    for dir in src.GetListOfKeys():
        print(dir.GetName())
        if (dir.GetName() == dir_name):
            continue
        cp_dir = tgt.mkdir(dir.GetName())
        src_dir = src.Get(dir.GetName())
        for obj in src_dir.GetListOfKeys():
            real_obj = src_dir.Get(obj.GetName())
            cp_dir.cd()
            real_obj.Write()  # Write each object to the target file

    # Remove the directory from the source file
    tgt.Delete(f"{dir_name};*")
    
    # Close files
    tgt.Close()
    src.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_file", "-i", metavar="text",
                        help="Input Root File", required=True)
    parser.add_argument("--output_file", "-o", metavar="text",
                        default="AnalysisResults_LHC23_Dstar_train293534.root",
                        help="output file", required=False)
    parser.add_argument("--target_dir", "-d", metavar="text",
                        default="hf-task-dplus",
                        help="output file", required=False)
    args = parser.parse_args()


    remove_directory(args.input_file, args.output_file, args.target_dir)