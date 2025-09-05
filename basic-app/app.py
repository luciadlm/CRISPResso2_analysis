import os
import re
import csv
import tempfile
import zipfile
import matplotlib.pyplot as plt

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
    ui.output_data_frame("results"),
    ui.output_plot("plot")
)

# Define Server
def server(input, output, session):

    @reactive.Calc
    def extracted_folder():
        if not input.zip_file():
            return None
        # Access the path of the uploaded zip file
        uploaded_zip = input.zip_file()[0]['datapath']
        print(uploaded_zip)
        tmpdir = tempfile.TemporaryDirectory()
        with zipfile.ZipFile(uploaded_zip, 'r') as zip_ref:
            zip_ref.extractall(tmpdir.name)
        return tmpdir  # return TemporaryDirectory object

    def process_treatment_folder(treatment_path, wt_seq, mut_seq):
        # Paths to expected CRISPResso files
        map_stat_path = os.path.join(treatment_path, "CRISPResso_mapping_statistics.txt")
        quant_edit_path = os.path.join(treatment_path, "CRISPResso_quantification_of_editing_frequency.txt")
        frameshift_path = os.path.join(treatment_path, "Frameshift_analysis.txt")
        amplicon_path = os.path.join(treatment_path, "CRISPResso_RUNNING_LOG.txt")
        alleles_path = os.path.join(treatment_path, "Alleles_frequency_table_around_sgRNA_AAAAAATGTTGGCCTCTCTT.txt")

        # Check if files exist
        if not all(os.path.exists(p) for p in [map_stat_path, quant_edit_path, frameshift_path, alleles_path]):
            return None

        # Use your functions to parse data
        total_reads = script.get_total_reads(map_stat_path)
        modif_data = script.get_modif_reads(quant_edit_path, total_reads)
        frame_data = script.get_frame_reads(frameshift_path, total_reads)
        amplicon_data = script.get_input_run(amplicon_path)
        mut_wt_data = script.get_mut_wt_reads(alleles_path, wt_seq, mut_seq, amplicon_data["amplicon_seq"], total_reads)

        return [
            os.path.basename(treatment_path),
            modif_data["wt_unmodified_percent"],
            modif_data["unmodified_reads"],
            mut_wt_data["mut_percent"],
            mut_wt_data["mut_reads"],
            mut_wt_data["wt_percent"],
            mut_wt_data["wt_reads"],
            frame_data["frameshift_percent"],
            frame_data["frameshift_reads"],
            frame_data["inframe_percent"],
            frame_data["inframe_reads"],
            modif_data["indel_percent"],
            modif_data["modified_reads"],
            mut_wt_data["mut_wt_percent"],
            total_reads,
        ]

    def get_replicate_data(tmpdir_path, wt_seq, mut_seq):
        results = []
        # print(tmpdir_path)
        for main_folder in os.listdir(tmpdir_path):
            main_folder_path = os.path.join(tmpdir_path, main_folder)
            # print(main_folder_path)
            if not os.path.isdir(main_folder_path):
                continue

            for replicate_name in os.listdir(main_folder_path):
                replicate_path = os.path.join(main_folder_path, replicate_name)
                # print(replicate_path)
                if not os.path.isdir(replicate_path):
                    continue
                
                denom_value = None
                
                for treatment_name in os.listdir(replicate_path):
                    treatment_path = os.path.join(replicate_path, treatment_name)
                    # print(treatment_path)
                    if not os.path.isdir(treatment_path):
                        continue

                    data_row = process_treatment_folder(treatment_path, wt_seq, mut_seq)
                    if data_row:
                        if denom_value is None and "Day02" in treatment_name:
                            denom_value = data_row[13]
                    
                        #if denom_value is None:
                            #denom_value = 1

                        sensitivity_value = round(data_row[13]/denom_value * 100, 2)
                        data_row.append(sensitivity_value)
                        data_row.insert(0, replicate_name)
                        results.append(data_row)

        return results

    @reactive.Calc
    def get_dataframe():
        import pandas as pd

        tmpdir = extracted_folder()
        #print(tmpdir.name)
        if tmpdir is None:
            return pd.DataFrame([["Please upload a zip file."]], columns=["Message"])

        mut_seq = input.mut_seq()
        wt_seq = input.wt_seq()
        if not mut_seq or not wt_seq:
            return pd.DataFrame([["Please enter both MUT and WT* sequences."]], columns=["Message"])

        # Write header
        header = ['Replicate', 'Sample set', 'WT/Unmodified %', 'Reads (WT/Unmodified)', 'MUT %', 'Reads (MUT)', 
                  'WT* %', 'Reads (WT*)', 'Frame shift %', 'Reads (Frame shift)', 'In-frame %', 'Reads (In-frame)', 
                  'Indel %', 'Reads (Modified)', 'MUT/WT* %', 'Total Reads', 'Sensitivity']

        data = get_replicate_data(tmpdir.name, wt_seq, mut_seq)

        if not data:
            return pd.DataFrame([["No valid results in zip file."]], columns=["Error"])
        
        df = pd.DataFrame(data, columns=header)
        return df 

    @output()
    @render.data_frame
    def results():
        return get_dataframe()
    
    @output()
    @render.plot
    def plot():
        df = get_dataframe()
        
        # Check if expected columns exist
        if "Sensitivity" not in df.columns or "Sample set" not in df.columns:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No data to plot", ha='center', va='center')
            ax.axis('off')
            return fig

        sensitivity = df["Sensitivity"]
        labels = df["Sample set"]

        fig, ax = plt.subplots(figsize=(10,6))
        ax.bar(labels, sensitivity)
        ax.set_title("Sensitivity by Sample set")
        ax.set_ylabel("Sensitivity")
        ax.set_xlabel("Sample set")

        return fig

app = App(app_ui, server)