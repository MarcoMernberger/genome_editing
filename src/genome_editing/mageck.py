import subprocess
from pathlib import Path
from typing import List


class Mageck:
    def __init__(self, output_folder: str = "results/crispresso"):
        self.output_folder = Path(output_folder)
        self._docker_command = [
            "docker",
            "run",
            "-v",
            "${ANYSNAKE2_PROJECT_DIR}:/project",
            "-w",
            "/project",
            "-i",
            "lucapinello/mageck",
        ]

    @property
    def docker_command(self):
        return self._docker_command

    @docker_command.setter
    def docker_command(self, docker_command: List[str]):
        self._docker_command = docker_command

    def call_mageck(
        self,
        # report_name: str,
        # input_files: List[str],
        # df_amplicons: Union[DataFrame, Callable],
        # output_folder: Path,
        # options: Optional[Dict[Any, Any]] = None,
        # quantification_window_size=10,
        # quantification_window_center=-3,
        # sample: Optional[str] = None,
    ):
        # if not isinstance(df_amplicons, DataFrame):
        #     try:
        #         df_amplicons = df_amplicons()
        #     except TypeError:
        #         raise TypeError("df_amplicons must be a DataFrame or a Callable.")
        # print(df_amplicons.Sample)
        # if sample is not None:
        #     df_amplicons = df_amplicons[df_amplicons.Sample == sample]
        #     if df_amplicons.shape[0] == 0:
        #         print(df_amplicons.Sample)
        #         print(sample)
        #         raise ValueError(f"Could not find {sample} in df_amplicons")
        # print(sample)
        # amplicon_seq = ",".join(df_amplicons.Amplicon.values)
        # sgrnas = ",".join(df_amplicons.sgRNA.values)
        # amplicon_names = ",".join(df_amplicons.Sample.values)
        command = self.docker_command + ["mageck", "-h"]
        # if len(input_files) == 2:
        #     command.extend(
        #         [
        #             "--fastq_r2",
        #             str(input_files[1]),
        #         ]
        #     )
        # if options is not None:
        #     command.append(dict_to_string_of_items(options))
        cmd = " ".join(command)
        print(cmd)
        try:
            subprocess.run(cmd, shell=True)
        except subprocess.CalledProcessError:
            print(cmd)
            raise


#     def call_crispresso2_PE(
#         self,
#         report_name: str,
#         input_files: List[str],
#         df_amplicons: Union[DataFrame, Callable],
#         output_folder: Path,
#         options: Optional[Dict[Any, Any]] = None,
#         sample: Optional[str] = None,
#         # quantification_window_size=10,
#         # quantification_window_center=-3,
#     ):
#         if not isinstance(df_amplicons, DataFrame):
#             try:
#                 df_amplicons = df_amplicons()
#             except TypeError:
#                 raise TypeError("df_amplicons must be a DataFrame or a Callable.")
#         if sample is not None:
#             df_amplicons = df_amplicons[df_amplicons.Sample == sample]
#             if df_amplicons.shape[0] == 0:
#                 print(df_amplicons.Sample)
#                 print(sample)
#                 raise ValueError(f"Could not find {sample} in df_amplicons")
#         print(sample)
#         print("df_amplicons")
#         print(df_amplicons)
#         amplicon_seq = ",".join(df_amplicons.Amplicon.values)
#         amplicon_names = ",".join(df_amplicons.name.values)
#         prime_editing_pegRNA_spacer_seq = ",".join(df_amplicons.pegRNAspacer.values)
#         prime_editing_pegRNA_extension_seq = ",".join(
#             df_amplicons.pegRNAextension.values
#         )
#         prime_editing_pegRNA_scaffold_seq = ",".join(df_amplicons.pegRNAscaffold.values)
#         # nick_guide = ",".join(df_amplicons.nick.values)
#         output_folder.mkdir(parents=True, exist_ok=True)
#         command = self.docker_command + [
#             "CRISPResso",  # , "-h"]
#             "--fastq_r1",
#             str(input_files[0]),
#         ]
#         if len(input_files) == 2:
#             command.extend(
#                 [
#                     "--fastq_r2",
#                     str(input_files[1]),
#                 ]
#             )
#         command += [
#             "--amplicon_seq",
#             amplicon_seq,
#             "--prime_editing_pegRNA_spacer_seq",
#             f"{prime_editing_pegRNA_spacer_seq}",
#             "--prime_editing_pegRNA_extension_seq",
#             f"{prime_editing_pegRNA_extension_seq}",
#             "--prime_editing_pegRNA_scaffold_seq",
#             f"{prime_editing_pegRNA_scaffold_seq}",
#             # "--prime_editing_nicking_guide_seq",
#             # f"{nick_guide}",
#             "--write_cleaned_report",
#             "--output_folder",
#             str(output_folder),
#             "--name",
#             report_name,
#             "--amplicon_name",
#             f"{amplicon_names}",
#             "--debug",
#         ]
#         # if raw_sample.is_paired:
#         if options is not None:
#             command.append(dict_to_string_of_items(options))
#         cmd = " ".join(command)
#         print(cmd)
#         try:
#             subprocess.run(cmd, shell=True)
#         except subprocess.CalledProcessError:
#             print(cmd)
#             raise

#     def run(
#         self,
#         report_name: str,
#         input_files: List[str],
#         df_amplicons: Union[DataFrame, Callable],
#         additional_folder: Optional[Path] = None,
#         options: Optional[Dict[Any, Any]] = None,
#         quantification_window_size: int = 10,
#         quantification_window_center: int = -3,
#         mode="classic",
#         dependencies: List[Job] = [],
#         sample: Optional[str] = None,
#     ):
#         output_folder = self.output_folder
#         if additional_folder is not None:
#             output_folder = output_folder / additional_folder
#         output_folder.mkdir(parents=True, exist_ok=True)
#         filename = f"CRISPResso_on_{report_name}.html"
#         outfile = output_folder / filename

#         def __dump(output_file):
#             if mode == "classic":
#                 self.call_crispresso2(
#                     report_name,
#                     input_files,
#                     df_amplicons,
#                     output_folder,
#                     options,
#                     quantification_window_size,
#                     quantification_window_center,
#                     sample=sample,
#                 )
#             elif mode == "prime_editing":
#                 self.call_crispresso2_PE(
#                     report_name,
#                     input_files,
#                     df_amplicons,
#                     output_folder,
#                     options,
#                     # quantification_window_size,
#                     # quantification_window_center,
#                     sample=sample,
#                 )
#             else:
#                 raise ValueError("mode must be 'classic' or 'prime_editing'.")

#         job = ppg.FileGeneratingJob(outfile, __dump).depends_on(dependencies)
#         job.cores_needed = 17  # we must do one after another because of docker ...
#         return job


# def call_crispresso2(
#     report_name: str,
#     input_files: List[str],
#     df_amplicons: DataFrame,
#     output_folder: Path,
#     options: Optional[Dict[Any, Any]] = None,
#     quantification_window_size=10,
#     quantification_window_center=-3,
# ):
#     amplicon_seq = ",".join(df_amplicons.Amplicon.values)
#     sgrnas = ",".join(df_amplicons.sgRNA.values)
#     amplicon_names = ",".join(df_amplicons.Gen.values)
#     output_folder.mkdir(parents=True, exist_ok=True)
#     docker_command = [
#         "docker",
#         "run",
#         "-v",
#         "${ANYSNAKE2_PROJECT_DIR}:/project",
#         "-w",
#         "/project",
#         "-i",
#         "pinellolab/crispresso2",
#     ]
#     command = docker_command + [
#         "CRISPResso",  # , "-h"]
#         "--fastq_r1",
#         str(input_files[0]),
#         "--amplicon_seq",
#         amplicon_seq,
#         "--guide_seq",
#         sgrnas,
#         "--quantification_window_size",
#         f"{quantification_window_size}",
#         "--quantification_window_center",
#         f"{quantification_window_center}",
#         "--base_editor_output",
#         "--exclude_bp_from_right",
#         "1",
#         "--exclude_bp_from_left",
#         "1",
#         "-o",
#         str(output_folder),
#         "-n",
#         report_name,
#         "-an",
#         amplicon_names,
#         "--debug",
#     ]
#     # if raw_sample.is_paired:
#     if len(input_files) == 2:
#         command.extend(
#             [
#                 "--fastq_r2",
#                 str(input_files[1]),
#             ]
#         )
#     if options is not None:
#         command.append(dict_to_string_of_items(options))
#     cmd = " ".join(command)
#     print(cmd)
#     try:
#         subprocess.run(cmd, shell=True)
#     except subprocess.CalledProcessError:
#         print(cmd)
#         raise
