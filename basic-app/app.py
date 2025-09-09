import os
import re
import csv
import tempfile
import zipfile
import glob
import matplotlib.pyplot as plt
import pandas as pd

from shiny import App, render, ui, reactive
#from shiny.express import input

import script

# Define UI
app_ui = ui.page_fluid(
    ui.panel_title("CRISPResso2 Analysis"),

    ui.input_file("zip_file", "Upload your zipped CRISPResso folders (.zip)", accept=[".zip"]),

    # Inputs: MUT and WT* sequences (text area)
    ui.input_text_area("mut_seq", "Enter the MUT sequence:", rows=3),
    ui.input_text_area("wt_seq", "Enter the WT* sequence:", rows=3),

    # Output: results text (TSV table as plain text)
    ui.layout_columns(
        ui.output_data_frame("results"),
        ui.output_plot("plot"),
        col_widths=[8, 4],
    ),
)

# Define Server
def server(input, output, session):

    @reactive.Calc
    def extracted_folder():
        """Extracts the uploaded ZIP file into a temporary directory."""

        if not input.zip_file():
            return None
        
        # Get the path of the uploaded zip file
        uploaded_zip = input.zip_file()[0]['datapath']

        # Create a temporary directory to extract contents
        tmpdir = tempfile.TemporaryDirectory()

        # Extract all contents of the zip file to the temporary directory
        with zipfile.ZipFile(uploaded_zip, 'r') as zip_ref:
            zip_ref.extractall(tmpdir.name)

        return tmpdir  # Return TemporaryDirectory object

    def process_treatment_folder(treatment_path, wt_seq, mut_seq):
        """Process a single CRISPResso output folder and extracts key editing metrics."""

        # Define the expected CRISPResso output files paths
        map_stat_path = os.path.join(treatment_path, "CRISPResso_mapping_statistics.txt")
        quant_edit_path = os.path.join(treatment_path, "CRISPResso_quantification_of_editing_frequency.txt")
        frameshift_path = os.path.join(treatment_path, "Frameshift_analysis.txt")
        amplicon_path = os.path.join(treatment_path, "CRISPResso_RUNNING_LOG.txt")

        # Look for the allele frequency table (filename may vary by sgRNA sequence)
        pattern = os.path.join(treatment_path, "Alleles_frequency_table_around_sgRNA_*.txt")
        matches = glob.glob(pattern)
        if not matches:
            return None # Allele file not found
        alleles_path = matches[0]

        # Verify that all required files exist before processing
        if not all(os.path.exists(p) for p in [map_stat_path, quant_edit_path, frameshift_path, alleles_path]):
            return None

        # Use the functions to extract the relevant data
        total_reads = script.get_total_reads(map_stat_path)
        modif_data = script.get_modif_reads(quant_edit_path, total_reads)
        frame_data = script.get_frame_reads(frameshift_path, total_reads)
        amplicon_data = script.get_input_run(amplicon_path)
        mut_wt_data = script.get_mut_wt_reads(
            alleles_path, 
            wt_seq, 
            mut_seq, 
            amplicon_data["amplicon_seq"], 
            total_reads
        )

        # Return all the processed information in a dictionary
        return {
            "Treatment": os.path.basename(treatment_path),
            "WT/Unmodified %": modif_data["wt_unmodified_percent"],
            "Reads (WT/Unmodified)": modif_data["unmodified_reads"],
            "MUT %": mut_wt_data["mut_percent"],
            "Reads (MUT)": mut_wt_data["mut_reads"],
            "WT* %": mut_wt_data["wt_percent"],
            "Reads (WT*)": mut_wt_data["wt_reads"],
            "Frame shift %": frame_data["frameshift_percent"],
            "Reads (Frame shift)": frame_data["frameshift_reads"],
            "In-frame %": frame_data["inframe_percent"],
            "Reads (In-frame)": frame_data["inframe_reads"],
            "Indel %": modif_data["indel_percent"],
            "Reads (Modified)": modif_data["modified_reads"],
            "MUT/WT* %": mut_wt_data["mut_wt_percent"],
            "Total Reads": total_reads,
        }

    def get_all_folders(base_path, wt_seq, mut_seq):
        """Walks through all subdirectories in the base path to find CRISPResso output folders,
        processes each valid folder and returns a DataFrame with all compiled results"""

        results = []
        
        denom_value = None

        # print(tmpdir_path)
        for root, dirs, files in os.walk(base_path):

            expected_files = [
                "CRISPResso_mapping_statistics.txt",
                "CRISPResso_quantification_of_editing_frequency.txt",
                "Frameshift_analysis.txt"
                # "Alleles_frequency_table_around_sgRNA_AAAAAATGTTGGCCTCTCTT.txt"
            ]

            if all(f in files for f in expected_files):
                row = process_treatment_folder(root, wt_seq, mut_seq)
                if row:
                    row["Replicate"] = os.path.basename(os.path.dirname(root))
                    #row["Treatment"] = os.path.basename(root)
                    results.append(row)

        if not results:
            return pd.DataFrame()

        df = pd.DataFrame(results)
        
        if denom_value is None:
            denom_value = df["MUT/WT* %"].iloc[0]  
        
        df["Sensitivity"] = round(df["MUT/WT* %"] / denom_value * 100, 2)

        return df


    @reactive.Calc
    def get_dataframe():
        
        tmpdir = extracted_folder()
        #print(tmpdir.name)
        if tmpdir is None:
            return pd.DataFrame([["Please upload a zip file."]], columns=["Message"])

        mut_seq = input.mut_seq()
        wt_seq = input.wt_seq()
        if not mut_seq or not wt_seq:
            return pd.DataFrame([["Please enter both MUT and WT* sequences."]], columns=["Message"])

        # Write header
        header = ['Replicate', 'Treatment', 'WT/Unmodified %', 'Reads (WT/Unmodified)', 'MUT %', 'Reads (MUT)', 
                  'WT* %', 'Reads (WT*)', 'Frame shift %', 'Reads (Frame shift)', 'In-frame %', 'Reads (In-frame)', 
                  'Indel %', 'Reads (Modified)', 'MUT/WT* %', 'Total Reads', 'Sensitivity']

        data = get_all_folders(tmpdir.name, wt_seq, mut_seq)

        if data.empty:
            return pd.DataFrame([["No valid results in zip file."]], columns=["Error"])
        
        df = pd.DataFrame(data, columns=header)
        return df 

    def sanitize_id(name):
        # Replace anything that is NOT a letter, number, or underscore with underscore
        return re.sub(r'[^a-zA-Z0-9_]', '_', name)

    @output()
    @render.data_frame
    def results():
        df = get_dataframe()
        return df
    
    @output()
    @render.plot
    def plot():
        df = get_dataframe()
        
        # Check if expected columns exist
        if "Sensitivity" not in df.columns or "Treatment" not in df.columns:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No data to plot", ha='center', va='center')
            ax.axis('off')
            return fig

        sensitivity = df["Sensitivity"]
        labels = df["Treatment"]

        fig, ax = plt.subplots(figsize=(10,6))
        ax.bar(labels, sensitivity)
        ax.set_title("Sensitivity by Treatment")
        ax.set_ylabel("Sensitivity")
        ax.set_xlabel("Treatment")

        return fig

app = App(app_ui, server)