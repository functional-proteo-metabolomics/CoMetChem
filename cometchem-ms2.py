import pastaq
import json
import os
import matplotlib
# NOTE: Using the default backend create problems when multithreaded with wxWidgets.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import wx
import threading
import sys

# Set the font size for the plots.
plt.rcParams.update({'font.size': 22})

def gauss(x, a, mu, sigma):
    return a * np.exp(-(x - mu)**2 / (2 * sigma**2))

def calculate_r2(x, y, fit):
    residuals = (y - gauss(x, *fit))
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res/ss_tot
    return r2

def theoretical_fwhm(mz, resolution):
    mz_ref = 200
    fwhm_ref = mz_ref/resolution
    fwhm = fwhm_ref * ((mz/mz_ref) ** 1.5)
    return fwhm

def calculate_light_heavy_ratio(x):
    if x.shape[0] != 2:
        return
    x = x.sort_values(by=['fragment_mz'])
    ratio_name = '{} / {}'.format(x['fragment_sequence'].iloc[0], x['fragment_sequence'].iloc[1])
    raw_height = x['raw_height'].iloc[0] / x['raw_height'].iloc[1]
    raw_total_intensity = x['raw_total_intensity'].iloc[0] / x['raw_total_intensity'].iloc[1]
    fitted_height = x['fitted_height'].iloc[0] / x['fitted_height'].iloc[1]
    fitted_area = x['fitted_area'].iloc[0] / x['fitted_area'].iloc[1]
    ratios = pd.Series({
        'ratio_name': ratio_name,
        'raw_height': raw_height,
        'raw_total_intensity': raw_total_intensity,
        'fitted_height': fitted_height,
        'fitted_area': fitted_area,
        })
    return ratios

def search_fragments(scan, target_fragments, out_dir, pastaq_parameters, plot):
    intensities = np.array(scan.intensity)
    mzs = np.array(scan.mz)

    # Create output directory if it doesn't exist.
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Prepare the results data frame.
    results = target_fragments.copy()
    results = results.reset_index(drop=True)
    results['raw_height'] = np.nan
    results['raw_total_intensity'] = np.nan
    results['fitted_height'] = np.nan
    results['fitted_mz'] = np.nan
    results['fitted_sigma_mz'] = np.nan
    results['fitted_area'] = np.nan
    results['fitted_r2'] = np.nan

    # Create a stemplot of the full MSMS scan.
    if plot:
        scan_fig, scan_ax = plt.subplots()
        scan_ax.ticklabel_format(useOffset=False)
        scan_ax.stem(
                mzs, intensities, use_line_collection=True, basefmt=' ', markerfmt=' ')

    # Search and quantify fragments from the inclusion list on the MSMS scan.
    for i in range(0, target_fragments.shape[0]):
        fragment_mz = target_fragments['fragment_mz'].iloc[i]
        fragment_ion = target_fragments['fragment_ion'].iloc[i]
        fragment_sequence = target_fragments['fragment_sequence'].iloc[i]

        # Mark the full MSMS scan stemplot with diamonds on the search sites.
        if plot:
            scan_ax.scatter(
                    x=fragment_mz,
                    y=0,
                    marker='d', s=30, color='crimson', zorder=1e5)

        # Check if there is data in the scan for this fragment.
        fragment_theoretical_sigma = theoretical_fwhm(
                fragment_mz, pastaq_parameters['resolution_msn']) / 2.355
        min_mz = fragment_mz - fragment_theoretical_sigma * 3
        max_mz = fragment_mz + fragment_theoretical_sigma * 3
        idx = (mzs > min_mz) & (mzs < max_mz)
        if idx.sum() == 0:
            continue

        # Perform Gaussian curve fitting for this fragment.
        x = mzs[idx]
        y = intensities[idx]
        a0 = np.max(y)
        mu0 = fragment_mz
        sigma0 = fragment_theoretical_sigma
        weights = 1 / gauss(x, 1, mu0, sigma0)
        try:
            fit, cov = curve_fit(
                        gauss, x, y,
                        p0=[a0, mu0, sigma0],
                        sigma=weights,
                        bounds=(
                            [0, fragment_mz - fragment_theoretical_sigma/2, fragment_theoretical_sigma / 2],
                            [np.inf, fragment_mz + fragment_theoretical_sigma/2, fragment_theoretical_sigma * 2],
                            )
                        )
        except Exception as e:
            print(e)
            continue
        a_fit, mu_fit, sigma_fit = fit
        x_fit = np.linspace(min_mz, max_mz, 100)
        y_fit = gauss(x_fit, *fit)

        # Save the quantification to the results table.
        raw_height = a0
        raw_total_intensity = y.sum()
        fitted_area = np.trapz(y)
        fitted_r2 = calculate_r2(x, y, fit)
        if fitted_r2 < pastaq_parameters['min_r2']:
            continue

        results.at[i, 'raw_height'] = raw_height
        results.at[i, 'raw_total_intensity'] = raw_total_intensity
        results.at[i, 'fitted_height'] = a_fit
        results.at[i, 'fitted_mz'] = mu_fit
        results.at[i, 'fitted_sigma_mz'] = sigma_fit
        results.at[i, 'fitted_area'] = fitted_area
        results.at[i, 'fitted_r2'] = fitted_r2

        # Plot the scan and fitting.
        if plot:
            fragment_fig, fragment_ax = plt.subplots()
            fragment_ax.ticklabel_format(useOffset=False)
            fragment_ax.stem(
                    x, y, use_line_collection=True, basefmt=' ', markerfmt=' ')
            fragment_ax.plot(x_fit, y_fit, linestyle='--', label='$R^2$ = {:.4f}'.format(fitted_r2))
            fragment_ax.set_xlabel('m/z')
            fragment_ax.set_ylabel('Intensity')
            fragment_ax.set_title('{}  {}'.format(fragment_ion, fragment_sequence))
            fragment_ax.legend()
            fragment_fig.set_size_inches(13 * 16/9, 13)
            fragment_fig.tight_layout()
            fragment_plot_name = os.path.join(
                    # out_dir, '{:03d}_{}_{:.5f}.svg'.format(i, fragment_ion, fragment_mz))
                    out_dir, '{:03d}_{}_{:.5f}.png'.format(i, fragment_ion, fragment_mz))
            plt.figure(fragment_fig.number)
            plt.savefig(fragment_plot_name, dpi=300)
            plt.close(fragment_fig)

    # Save the results of the fragment search to csv.
    results_file_name = os.path.join(out_dir, 'fragments.csv')
    results.to_csv(results_file_name, index=False)

    # Create ratios between light/heavy ions.
    light_heavy_ratio = results.groupby('fragment_ion').apply(calculate_light_heavy_ratio)
    light_heavy_ratio_file_name = os.path.join(out_dir, 'light_heavy_ratio.csv')
    light_heavy_ratio.to_csv(light_heavy_ratio_file_name, index=True)

    # Save the scan stemplot.
    if plot:
        scan_ax.set_xlabel('m/z')
        scan_ax.set_ylabel('Intensity')
        scan_fig.set_size_inches(13 * 16/9, 13)
        scan_fig.tight_layout()
        scan_plot_name = os.path.join(out_dir, 'scan.svg')
        plt.figure(scan_fig.number)
        plt.savefig(scan_plot_name, dpi=300)
        plt.close(scan_fig)

def run_analysis(i, sample_name, target_samples, target_fragments, mzml_dir, output_dir, pastaq_parameters, create_plots, override_existing):
    file_name = '{}.mzML'.format(os.path.join(mzml_dir, sample_name))
    if not os.path.exists(file_name):
        print("Can't open file:", file_name)
        return

    # Read raw data and store it as binary for future use.
    out_path = os.path.join(output_dir, 'raw', "{}.ms2".format(sample_name))
    if not os.path.exists(out_path) or override_existing:
        print("Reading MS2 mzML raw data for:", sample_name)
        raw_data = pastaq.read_mzml(
                file_name,
                min_mz=pastaq_parameters['min_mz'],
                max_mz=pastaq_parameters['max_mz'],
                min_rt=pastaq_parameters['min_rt'],
                max_rt=pastaq_parameters['max_rt'],
                instrument_type=pastaq_parameters['instrument_type'],
                resolution_ms1=pastaq_parameters['resolution_ms1'],
                resolution_msn=pastaq_parameters['resolution_msn'],
                reference_mz=pastaq_parameters['reference_mz'],
                fwhm_rt=pastaq_parameters['avg_fwhm_rt'],
                polarity=pastaq_parameters['polarity'],
                ms_level=2,
                )
        print("Saving MS2 raw data for:", sample_name, "-", out_path)
        raw_data.dump(out_path)
    else:
        print("Reading MS2 raw data for:", sample_name)
        raw_data = pastaq.read_raw_data(out_path)

    if 'precursor_min_rt' in target_samples.columns:
        min_rt = target_samples['precursor_min_rt'].iloc[i] * 60.0
    else:
        precursor_rt = target_samples['precursor_rt'].iloc[i] * 60.0
        min_rt = precursor_rt - 1.5 * raw_data.fwhm_rt
    if 'precursor_max_rt' in target_samples.columns:
        max_rt = target_samples['precursor_max_rt'].iloc[i] * 60.0
    else:
        precursor_rt = target_samples['precursor_rt'].iloc[i] * 60.0
        max_rt = precursor_rt + 1.5 * raw_data.fwhm_rt

    precursor_mz = target_samples['precursor_mz'].iloc[i]
    min_mz = precursor_mz - 1.5 * theoretical_fwhm(precursor_mz, pastaq_parameters['resolution_ms1'])
    max_mz = precursor_mz + 1.5 * theoretical_fwhm(precursor_mz, pastaq_parameters['resolution_ms1'])
    print('Finding precursors within min_rt: {:.2f}, max_rt: {:.2f}, min_mz: {:.5f}, max_mz: {:.5f}'.format(min_rt, max_rt, min_mz, max_mz))
    num_found = 0
    for scan in raw_data.scans:
        if (scan.retention_time > min_rt and
                scan.retention_time < max_rt and
                scan.precursor_information.mz > min_mz and
                scan.precursor_information.mz < max_mz):
            num_found += 1
            fragment_table = target_fragments[
                    (target_fragments['precursor_mz'] > min_mz) &
                    (target_fragments['precursor_mz'] < max_mz)]

            out_dir_results = os.path.join(output_dir, 'results', '{}_{:.5f}_{:.2f}__{:.2f}_{}'.format(
                sample_name, precursor_mz, min_rt, max_rt, scan.scan_number))
            search_fragments(scan, fragment_table, out_dir_results, pastaq_parameters, create_plots)
    print("Found:", num_found)

class MainWindow(wx.Frame):
    def __init__(self):
        # Initialize window and main sizer.
        wx.Frame.__init__(self, None, title='MS2 Quantifier')
        self.panel = wx.Panel(self)
        global_sizer_v = wx.BoxSizer(wx.VERTICAL)

        # Initialize default parameters.
        self.generate_plots = True
        self.override_existing = False
        self.pastaq_parameters = pastaq.default_parameters('orbitrap', 10.0)
        self.pastaq_parameters['resolution_ms1'] = 60000
        self.pastaq_parameters['resolution_msn'] = 7500
        self.pastaq_parameters['min_r2'] = 0.25
        self.output_dir = ''
        self.mzml_dir = ''
        self.target_samples = ''
        self.target_fragments = ''

        # Prepare the inclusion list file selection box.
        inc_list_box = wx.StaticBox(
                self.panel, -1, 'Please select the inclusion list tables',
                style=wx.ALIGN_CENTRE_HORIZONTAL)
        inc_list_sizer_v = wx.StaticBoxSizer(inc_list_box, wx.VERTICAL)
        inc_list_sizer_h = wx.BoxSizer(wx.HORIZONTAL)
        inc_list_fp_wildcard = "CSV files (*.csv)|*.csv"

        # Sample list file picker.
        self.label_1 = wx.StaticText(
                self.panel, style=wx.ALIGN_RIGHT, label="Sample list:")
        self.inc_list_fp_1 = wx.FilePickerCtrl(
                self.panel, message="Choose a file",
                style=wx.FLP_USE_TEXTCTRL | wx.FD_FILE_MUST_EXIST,
                wildcard=inc_list_fp_wildcard)
        self.inc_list_fp_1.Bind(wx.EVT_FILEPICKER_CHANGED, self.inc_list_on_file_open)
        inc_list_sizer_h.Add(self.label_1, 1, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        inc_list_sizer_h.Add(self.inc_list_fp_1, 5, wx.ALL | wx.CENTER | wx.EXPAND, 5)

        # Fragment list file picker.
        self.label_2 = wx.StaticText(
                self.panel, style=wx.ALIGN_RIGHT, label="Fragment list:")
        self.inc_list_fp_2 = wx.FilePickerCtrl(
                self.panel, message="Choose a file",
                style=wx.FLP_USE_TEXTCTRL | wx.FD_FILE_MUST_EXIST,
                wildcard=inc_list_fp_wildcard)
        self.inc_list_fp_2.Bind(wx.EVT_FILEPICKER_CHANGED, self.inc_list_on_file_open)
        inc_list_sizer_h.Add(self.label_2, 1, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        inc_list_sizer_h.Add(self.inc_list_fp_2, 5, wx.ALL | wx.CENTER | wx.EXPAND, 5)

        inc_list_sizer_v.Add(inc_list_sizer_h, 1, wx.ALL | wx.CENTER | wx.EXPAND, 5)

        # Prepare the parameters sizer box.
        parameters_box = wx.StaticBox(
                self.panel, -1, 'Parameters',
                style=wx.ALIGN_CENTRE_HORIZONTAL)
        parameters_sizer_v = wx.StaticBoxSizer(parameters_box, wx.VERTICAL)
        parameters_sizer_h_1 = wx.BoxSizer(wx.HORIZONTAL)
        parameters_sizer_h_2 = wx.BoxSizer(wx.HORIZONTAL)
        parameters_sizer_h_checkboxes = wx.BoxSizer(wx.HORIZONTAL)

        # Resolution selection (MS1).
        self.resolution_label_1 = wx.StaticText(
                self.panel, style=wx.ALIGN_RIGHT, label="Resolution MS1:")
        self.resolution_num_1 = wx.TextCtrl(
                self.panel, -1, value=str(self.pastaq_parameters['resolution_ms1']))

        self.resolution_num_1.Bind(wx.EVT_CHAR, self.skip_non_numbers)
        self.resolution_num_1.Bind(wx.EVT_TEXT, self.update_ms1_resolution)
        parameters_sizer_h_1.Add(self.resolution_label_1, 1, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        parameters_sizer_h_1.Add(self.resolution_num_1, 5, wx.ALL | wx.CENTER | wx.EXPAND, 5)

        # Resolution selection (MS2).
        self.resolution_label_2 = wx.StaticText(
                self.panel, style=wx.ALIGN_RIGHT, label="Resolution MS2:")
        self.resolution_num_2 = wx.TextCtrl(
                self.panel, -1, value=str(self.pastaq_parameters['resolution_msn']))
        self.resolution_num_2.Bind(wx.EVT_CHAR, self.skip_non_numbers)
        self.resolution_num_2.Bind(wx.EVT_TEXT, self.update_msn_resolution)
        parameters_sizer_h_1.Add(self.resolution_label_2, 1, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        parameters_sizer_h_1.Add(self.resolution_num_2, 5, wx.ALL | wx.CENTER | wx.EXPAND, 5)

        # Input dir picker.
        self.input_dir_label = wx.StaticText(
                self.panel, style=wx.ALIGN_RIGHT, label="Input dir (mzML):")
        self.input_dir_fp = wx.DirPickerCtrl(self.panel, -1)
        self.input_dir_fp.Bind(wx.EVT_DIRPICKER_CHANGED, self.update_input_dir)
        parameters_sizer_h_2.Add(self.input_dir_label, 1, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        parameters_sizer_h_2.Add(self.input_dir_fp, 5, wx.ALL | wx.CENTER | wx.EXPAND, 5)

        # Output dir picker.
        self.output_dir_label = wx.StaticText(
                self.panel, style=wx.ALIGN_RIGHT, label="Output dir:")
        self.output_dir_fp = wx.DirPickerCtrl(self.panel, -1)
        self.output_dir_fp.Bind(wx.EVT_DIRPICKER_CHANGED, self.update_output_dir)
        parameters_sizer_h_2.Add(self.output_dir_label, 1, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        parameters_sizer_h_2.Add(self.output_dir_fp, 5, wx.ALL | wx.CENTER | wx.EXPAND, 5)

        # Add checkboxes.
        self.generate_plots_cb = wx.CheckBox(self.panel, label="Generate plots", size=(-1,50))
        self.generate_plots_cb.SetValue(True)
        self.generate_plots_cb.Bind(wx.EVT_CHECKBOX, self.generate_plots_cb_click)
        parameters_sizer_h_checkboxes.Add(self.generate_plots_cb, 1, wx.ALL | wx.CENTER, 5)

        # Bundle parameters in vertical sizer.
        parameters_sizer_v.Add(parameters_sizer_h_1, 0, wx.ALL | wx.EXPAND , 5)
        parameters_sizer_v.Add(parameters_sizer_h_2, 0, wx.ALL | wx.EXPAND , 5)
        parameters_sizer_v.Add(wx.StaticLine(self.panel, -1), 0, wx.ALL | wx.EXPAND, 5)
        parameters_sizer_v.Add(parameters_sizer_h_checkboxes, 0, wx.ALL, 5)

        # Buttons.
        buttons_sizer_h = wx.BoxSizer(wx.HORIZONTAL)
        self.start_button = wx.Button(self.panel, -1, label="Start", size=(200,50))
        self.start_button.Bind(wx.EVT_BUTTON, self.start_button_click)
        self.load_parameters_button = wx.Button(self.panel, -1, label="Load parameters", size=(200,50))
        self.load_parameters_button.Bind(wx.EVT_BUTTON, self.load_parameters)
        self.save_parameters_button = wx.Button(self.panel, -1, label="Save parameters", size=(200,50))
        self.save_parameters_button.Bind(wx.EVT_BUTTON, self.save_parameters)
        buttons_sizer_h.Add(self.start_button, 1, wx.ALL | wx.CENTER, 5)
        buttons_sizer_h.Add(self.load_parameters_button, 1, wx.ALL | wx.CENTER, 5)
        buttons_sizer_h.Add(self.save_parameters_button, 1, wx.ALL | wx.CENTER, 5)

        # Configure the final layout.
        global_sizer_v.Add(inc_list_sizer_v, 0, wx.ALL | wx.CENTER | wx.EXPAND, 5)
        global_sizer_v.Add(parameters_sizer_v, 0, wx.ALL | wx.CENTER | wx.EXPAND, 5)
        global_sizer_v.Add(buttons_sizer_h, 0, wx.ALL | wx.CENTER, 5)
        self.panel.SetSizer(global_sizer_v)
        self.Center()
        self.panel.Fit()
        self.Show()

    def inc_list_on_file_open(self, event):
        self.target_samples = self.inc_list_fp_1.GetPath()
        self.target_fragments = self.inc_list_fp_2.GetPath()

    def generate_plots_cb_click(self, event):
        self.generate_plots = not self.generate_plots

    def skip_non_numbers(self, event):
        keycode = event.GetKeyCode()
        key = chr(keycode)
        if key.isnumeric() or key in (
                'KEY_BACKSPACE', '\b', '\x7f', '\t',
                chr(wx.WXK_LEFT), chr(wx.WXK_RIGHT),
                chr(wx.WXK_CONTROL_C), chr(wx.WXK_CONTROL_V), chr(wx.WXK_CONTROL_X),
                ):
            event.Skip()

    def update_ms1_resolution(self, event):
        value = self.resolution_num_1.GetValue()
        if not value:
            self.pastaq_parameters['resolution_ms1'] = 0
        elif value.isnumeric():
            self.pastaq_parameters['resolution_ms1'] = int(self.resolution_num_1.GetValue())
        else:
            self.resolution_num_1.SetValue(str(self.pastaq_parameters['resolution_ms1']))

    def update_msn_resolution(self, event):
        value = self.resolution_num_2.GetValue()
        if not value:
            self.pastaq_parameters['resolution_msn'] = 0
        elif value.isnumeric():
            self.pastaq_parameters['resolution_msn'] = int(self.resolution_num_2.GetValue())
        else:
            self.resolution_num_2.SetValue(str(self.pastaq_parameters['resolution_msn']))

    def update_input_dir(self, event):
        self.mzml_dir = event.GetPath()

    def update_output_dir(self, event):
        self.output_dir = event.GetPath()

    def save_parameters(self, event):
        params = self.pastaq_parameters
        params['output_dir'] = self.output_dir
        params['mzml_dir'] = self.mzml_dir
        params['generate_plots'] = self.generate_plots
        params['override_existing'] = self.override_existing
        params['target_samples'] = self.target_samples
        params['target_fragments'] = self.target_fragments

        dlg = wx.DirDialog(self, "Choose a directory:", style=wx.DD_DEFAULT_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            file_name = os.path.join(dlg.GetPath(), 'parameters.json')
            with open(file_name, 'w') as json_file:
                json.dump(params, json_file)
            msg_dlg = wx.MessageDialog(
                    self, 'Parameters saved at:\n{}'.format(file_name), 'Success',
                    wx.OK|wx.STAY_ON_TOP)
            msg_dlg.ShowModal()
            msg_dlg.Destroy()
        dlg.Destroy()

    def load_parameters(self, event):
        wildcard = "Parameter files (*.json)|*.json"
        dlg = wx.FileDialog(
                self, "Choose the parameters file:",
                style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST, wildcard=wildcard)
        if dlg.ShowModal() == wx.ID_OK:
            file_name = dlg.GetPath()

        with open(file_name) as json_file:
            params = json.load(json_file)

            # Load parameters into internal data structures.
            self.pastaq_parameters = params
            self.output_dir = params['output_dir']
            self.mzml_dir = params['mzml_dir']
            self.generate_plots = params['generate_plots']
            self.override_existing = params['override_existing']
            self.target_samples = params['target_samples']
            self.target_fragments = params['target_fragments']

            # Update the GUI.
            self.inc_list_fp_1.SetPath(self.target_samples)
            self.inc_list_fp_2.SetPath(self.target_fragments)
            self.resolution_num_1.SetValue(str(self.pastaq_parameters['resolution_ms1']))
            self.resolution_num_2.SetValue(str(self.pastaq_parameters['resolution_msn']))
            self.input_dir_fp.SetPath(self.mzml_dir)
            self.output_dir_fp.SetPath(self.output_dir)
            self.generate_plots_cb.SetValue(self.generate_plots)

            msg_dlg = wx.MessageDialog(
                    self, 'Loaded parameters from:\n{}'.format(file_name), 'Success',
                    wx.OK|wx.STAY_ON_TOP)
            msg_dlg.ShowModal()
            msg_dlg.Destroy()

    def update_progress(self, i, msg):
        if self.running:
            (self.running, skip) = self.progress_dlg.Update(i, msg)

    def handle_progressbar_click(self, event):
        self.running = False
        self.work_thread.join()
        self.progress_dlg.Close()
        self.progress_dlg.Destroy()

    def run(self, target_samples, target_fragments):
        for i, sample_name in enumerate(target_samples['sample_name']):
            if not self.running:
                break
            msg = '{} ({} m/z)'.format(
                    sample_name, target_samples['precursor_mz'].iloc[i])
            run_analysis(
                i, sample_name,
                target_samples, target_fragments, self.mzml_dir,
                self.output_dir, self.pastaq_parameters, self.generate_plots,
                self.override_existing)
            wx.CallAfter(self.update_progress, i, msg)
            # FIXME: Need to wait 100ms because unknown race condition causes segfault.
            wx.MilliSleep(100)

        # FIXME: Need to wait 100ms because unknown race condition causes segfault.
        wx.MilliSleep(100)
        wx.CallAfter(self.update_progress, target_samples.shape[0], 'Done!')

    def start_button_click(self, event):
        # Read inclusion list data.
        target_samples = pd.read_csv(self.target_samples)
        target_fragments = pd.read_csv(self.target_fragments)
        output_dir = self.output_dir

        # Crete the output directories and subdirectories.
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if not os.path.exists(os.path.join(output_dir, 'raw')):
            os.makedirs(os.path.join(output_dir, 'raw'))

        max_samples = target_samples.shape[0]
        self.progress_dlg = wx.ProgressDialog(
            'Analysis progress', 'Analysis progress', maximum=max_samples, parent=self,
            style = wx.PD_CAN_ABORT|wx.PD_APP_MODAL)
        self.progress_dlg.Bind(wx.EVT_BUTTON, self.handle_progressbar_click)
        self.running = True

        self.work_thread = threading.Thread(
                target=self.run, args=[target_samples, target_fragments])

        self.work_thread.setDaemon(True)
        self.work_thread.start()

if __name__ == "__main__":
    app = wx.App()
    frame = MainWindow()
    app.MainLoop()
    sys.exit(0)
