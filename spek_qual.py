#!/usr/bin/env python3
from pathlib import Path
import os
import argparse
import json
import warnings

import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.optimize import minimize_scalar


# TODO: better path handling needed
PATH2DAT = "data"
PATH2RES = "result"


def search_string_in_file(path, string_to_search):
	"""
	Search for the given string in file and return lines containing that string,
	along with line numbers

	:param path: path to the file
	:type path: str
	:param string_to_search: string to search
	:type string_to_search: str
	:return: list with string and number line
	:rtype: list
	"""
	line_number = 0
	list_of_results = []
	# Open the file in read only mode
	with open(path, 'r') as read_obj:
		# Read all lines in the file one by one
		for line in read_obj:
			# For each line, check if line contains the string
			line_number += 1
			if string_to_search in line:
				# If yes, then add the line number & line as a tuple in the list
				list_of_results.append((line_number, line.rstrip()))
	# Return list of tuples containing line numbers and lines where string is found
	return list_of_results


def dir_exists(path: str, mkdir=False, noprint=False):
	"""
	Checks if directory under path exists, optionally creates directories in path

	:param path: path to the directory
	:param mkdir: True to create directory (directories tree) under path
	:param noprint: print information while executing or not
	:type noprint: flag to print info
	:return: bool
	"""
	if os.path.isdir(path):
		if not noprint:
			print(f'Directory "{path}" exists...')
		return True
	else:
		if not noprint:
			print(f'Directory "{path}" does not exist!')
		if mkdir:
			if not noprint:
				print(f'Making "{path}" directory!')
			os.makedirs(path)
		return False


def file_exists(path, noprint=False):
	"""
	Checks if file under path exists

	:param path: path to the file
	:type path: str
	:param noprint: print information while executing or not
	:type noprint: flag to print
	:return: bool
	"""
	if os.path.isfile(path):
		if not noprint:
			print(f'File "{path}" exists...')
		return True
	else:
		if not noprint:
			print(f'File "{path}" does not exist!')
		return False


def load_json(path):
	"""
	Load data file in json format from path

	:return: loaded json file
	:rtype: dict
	"""
	with open(path, "r") as fh:
		return json.load(fh)


def split_fluka_data(
		path, dat_type='F', det_n=None, save=True, noprint=False):
	"""
	Splits FLUKA data on detectors to different files

	:param path: full path to data with FLUKA data
	:type path: str
	:param dat_type: type of data 'F' for fluence, 'Y' for yield
	:type dat_type: str
	:param det_n: number of detector in file to load (Detector n). By default,
		is None, which means every index
	:type det_n: int
	:param save: save file to "result" folder or not
	:type save: bool
	:param noprint: print information while executing or not
	:type noprint: bool
	:return: dictionary with DataFrames of detectors {"detname": DataFrame, ...}
	"""
	if dir_exists(PATH2DAT, noprint=noprint):  # Check if directory exists
		if file_exists(path, noprint=noprint):  # Check if file exists
			pass
		else:
			exit()
	else:
		exit()

	dir_exists(PATH2RES, mkdir=True, noprint=noprint)

	if dat_type == 'Y':
		names = [
			'Emin [GeV]', 'Emax [GeV]', 'Y [GeV-1 primary-1]', 'Y_rerr [%]']
	elif dat_type == "LET":
		names = [
			'LETmin [keV/um]', 'LETmax [keV/um]', 'Yield [keV/um-1 primary-1]',
			'F_rerr [%]']
	else:
		names = [
			'Emin [GeV]', 'Emax [GeV]', 'F [cm-2 GeV-1 primary-1]',
			'F_rerr [%]']

	detectors = {}  # Dict with detectors' DataFrames from file
	parse = search_string_in_file(path, '#')  # Parse all lines with "#"
	detn = len(
		search_string_in_file(path, 'Detector n:'))  # Number of detectors

	if det_n is None:  # Split all detectors by default
		for i in range(detn):
			line_n = parse[i * 2 + 1][0]  # how many rows to skip = line number
			detname = parse[i * 2][1].split()[4]
			# pprint(detname)
			nrows = int(parse[i * 2 + 1][1].split()[
							5])  # number of energy intervals = rows

			data = pd.read_csv(
				path, sep='\s+', skiprows=line_n, nrows=nrows,
				names=names, usecols=names, header=None)

			if save:
				name = path.split('/')[-1].split('.')[0]
				data.to_csv(
					f'result/{name}_{detname}.tsv',
					sep='\t', index=False, header=True)
			detectors[detname] = data  # TODO !!!
	else:
		i = det_n - 1
		if 1 <= det_n < detn + 1:
			line_n = parse[i * 2 + 1][0]  # how many rows to skip = line number
			detname = parse[i * 2][1].split()[4]
			nrows = int(parse[i * 2 + 1][1].split()[
							5])  # number of energy intervals = rows

			data = pd.read_csv(
				path, sep='\s+', skiprows=line_n, nrows=nrows,
				names=names, usecols=names, header=None)

			if save:
				name = path.split('/')[-1].split('.')[0]
				data.to_csv(
					f'result/{name}_{detname}.tsv',
					sep='\t', index=False, header=True)
			detectors[detname] = data  # TODO !!!
	return detectors


def load_fluka_data(filename, det_n, row_drop=2, save=True, noprint=True):
	"""
	Load FLUKA-like data as DataFrame
	Load output FLUKA data in *tab.lis format:

	| Emin [GeV] | Emax [GeV] | F [cm-2 GeV-1 primary-1] | F_rerr [%] |

	:param filename: name of the file to load spectrum from
	:type filename: str
	:param det_n: number of detector in file to load (Detector n)
	:type det_n: int
	:param row_drop: how many rows to drop from beginning (cut to 1 keV)
	:type row_drop: int
	:param save: save file to "result" folder or not
	:type save: bool
	:param noprint: print information while executing or not
	:type noprint: bool
	:return: DataFrame with loaded from FLUKA output file data
	:rtype: pd.DataFrame | None
	"""
	spectrum = split_fluka_data(
		f"data/spectra/{filename}", det_n=det_n, save=save, noprint=noprint)
	spectrum = spectrum[list(spectrum.keys())[0]]

	# Get rid of first rows with
	spectrum = spectrum.tail(spectrum.shape[0] - row_drop)

	spectrum["Emin [keV]"] = spectrum["Emin [GeV]"]*1e6
	spectrum["Emax [keV]"] = spectrum["Emax [GeV]"]*1e6
	spectrum["dE [keV]"] = spectrum["Emax [keV]"]-spectrum["Emin [keV]"]

	# Make normalized spectrum
	spectrum["S(E)"] = \
		spectrum["F [cm-2 GeV-1 primary-1]"]\
		/ spectrum["F [cm-2 GeV-1 primary-1]"].sum()
	spectrum["Emid [keV]"] = spectrum["Emin [keV]"]+spectrum["dE [keV]"]/2

	return spectrum[["Emin [keV]", "Emax [keV]", "Emid [keV]", "dE [keV]", "S(E)"]]


def load_spekpy_data(filename, sep='\t', char=False, noprint=True):
	"""
	Load spekpy-like data as DataFrame
	Load output FLUKA data in a format:

	| Emid [keV] | F_total [keV-1 cm-2] | F_char [keV-1 cm-2] |

	:param filename: name of the file to load spectrum from
	:type filename: str
	:param sep: delimiter used to separate columns in file
	:type sep: str
	:param char: flag to get column with characteristic fluence, otherwise
		total fluence will be considered
	:type char: bool
	:param noprint: print information while executing or not
	:type noprint: bool
	:return: DataFrame with loaded from file data
	:rtype: pd.DataFrame
	"""
	names = ['Emid [keV]', 'F_total [keV-1 cm-2]', 'F_char [keV-1 cm-2]']
	spectrum = pd.read_csv(
		f"data/spectra/{filename}", sep=sep, comment="#",
		names=names, usecols=names)

	if not noprint:
		print("Verify loaded data:")
		print(spectrum)

	# TODO: here we assume that bin width is constant!
	spectrum["dE [keV]"] = spectrum["Emid [keV]"][1]-spectrum["Emid [keV]"][0]
	spectrum["Emin [keV]"] = spectrum["Emid [keV]"]-spectrum["dE [keV]"]/2
	spectrum["Emax [keV]"] = spectrum["Emid [keV]"]+spectrum["dE [keV]"]/2

	# Make normalized spectrum, use total fluence or only characteristic
	if char:
		fluence_col = 'F_char [keV-1 cm-2]'
	else:
		fluence_col = "F_total [keV-1 cm-2]"
	spectrum["S(E)"] = spectrum[fluence_col]/spectrum[fluence_col].sum()

	return spectrum[["Emin [keV]", "Emax [keV]", "Emid [keV]", "dE [keV]", "S(E)"]]


def xvl(
		spectrum, att_material="z13", att_material_rho=2.699, red_fraction=None,
		mu_source="nist", noprint=True):
	"""
	A function to calculate the first, second HVL (half value layer),
	QVL and TVL (in mm) for the desired material (att_material) with density
	of att_material_rho

	:param spectrum: DataFrame in a format:
		| Emin [keV] | Emax [keV] | Emid [keV] | dE [keV] | S(E) |
	:type spectrum: pd.DataFrame
	:param att_material: name of the attenuation material from NIST database,
		could be also specified by fraction in dict {Z:frac, ...}
	:type att_material: str | dict
	:param att_material_rho: density of attenuation material
	:type att_material_rho: float
	:param red_fraction: fraction to return value for, if is None (by default)
		then return standard (HVL1, HVL2, QVL, TVL) values as tuple, otherwise
		only one layer for specified fraction
	:type red_fraction: None | float
	:param mu_source: Load mu data from "nist" or "penelope" databases
	:type mu_source: str
	:param noprint: print information while executing or not
	:type noprint: bool
	:return: result as tuple (HVL1, HVL2, QVL, TVL) if red_fraction is None
		or XVL for the specified fraction red_fraction = X
	:rtype: float
	"""
	def cost_function(thick):
		"""
		A cost function is a measure of how wrong the model is in terms of its
		ability to estimate the relationship between X and y.
		Here MSE(Mean Squared Error) is used as a metric of how close our value to
		the desired.

		:param thick: thickness of the attenuator is a free parameter
		:type thick: float
		:return: Mean Squared Error
		:rtype: float
		"""
		# get rid of overflow, invalid value warnings
		warnings.filterwarnings("ignore")
		spectrum["Exposure 2"] =\
			spectrum["Emid [keV]"] * \
			spectrum["S(E)"] * \
			spectrum["mu_en/rho (air) (cm^2/g)"] * \
			spectrum["dE [keV]"] * \
			np.exp(
				-spectrum[f"mu ({att_material_name}) (cm-1)"] * thick)
		exposure2 = spectrum["Exposure 2"].sum()
		mse = np.power(red_fraction - (exposure2 / exposure1), 2)
		if not noprint:
			print(mse)
		return mse

	# Get air mass energy-absorption coefficients
	mu_en_rho_air = GetMu(
		medium="air", rho=1.0, mu_type="energy", mu_source=mu_source)
	spectrum["mu_en/rho (air) (cm^2/g)"] = mu_en_rho_air(spectrum["Emid [keV]"])

	# Get linear coefficient of attenuation
	att_material_name = att_material if isinstance(att_material, str) else "mixed"
	mu_att = GetMu(
		medium=att_material, rho=att_material_rho,
		mu_type="mass", mu_source=mu_source)
	spectrum[f"mu ({att_material_name}) (cm-1)"] = mu_att(spectrum["Emid [keV]"])

	# Calculate initial exposure without attenuator
	spectrum["Exposure 1"] =\
		spectrum["Emid [keV]"] * \
		spectrum["S(E)"] * \
		spectrum["mu_en/rho (air) (cm^2/g)"] * \
		spectrum["dE [keV]"]
	exposure1 = spectrum["Exposure 1"].sum()

	bracket = (1e-1, 1)  # TODO: we need this constrain for monoenergies
	if red_fraction is None:
		# Calculate all the layers to reduce Kerma to 0.5, 0.25, 0.1
		# We need to minimize cost to get the most accurate result
		red_fraction = 0.5
		hvl1 = minimize_scalar(
			cost_function, bracket=bracket, method="brent").x*10  # to mm
		red_fraction = 0.25
		qvl = minimize_scalar(
			cost_function, bracket=bracket, method="brent").x*10  # to mm
		hvl2 = qvl - hvl1
		red_fraction = 0.1
		tvl = minimize_scalar(
			cost_function, bracket=bracket, method="brent").x*10  # to mm
		return hvl1, hvl2, qvl, tvl
	else:
		return minimize_scalar(
			cost_function, bracket=bracket, method="brent").x*10  # to mm


def attenuate(
		spectrum, thick, att_material="z13", att_material_rho=2.699,
		mu_source="nist", noprint=True):
	"""
	A function to attenuate spectrum on a specified material thickness

	:param spectrum: DataFrame in a format:
		| Emin [keV] | Emax [keV] | Emid [keV] | dE [keV] | S(E) |
	:type spectrum: pd.DataFrame
	:param thick: thickness of attenuation material in mm
	:type thick: float
	:param att_material: name of the attenuation material from NIST database,
		could be also specified by fraction in dict {Z:frac, ...}
	:type att_material: str | dict
	:param att_material_rho: density of attenuation material
	:type att_material_rho: float
	:param mu_source: Load mu data from "nist" or "penelope" databases
	:type mu_source: str
	:param noprint: print information while executing or not
	:type noprint: bool
	:return: result as fraction on specified thickness
	:rtype: float
	"""
	# Get air mass energy-absorption coefficients
	mu_en_rho_air = GetMu(
		medium="air", rho=1.0, mu_type="energy", mu_source=mu_source)
	spectrum["mu_en/rho (air) (cm^2/g)"] = mu_en_rho_air(spectrum["Emid [keV]"])

	# Get linear coefficient of attenuation
	att_material_name = att_material if isinstance(att_material, str) else "mixed"
	mu_att = GetMu(
		medium=att_material, rho=att_material_rho,
		mu_type="mass", mu_source=mu_source)
	spectrum[f"mu ({att_material_name}) (cm-1)"] = mu_att(spectrum["Emid [keV]"])

	if not noprint:
		print("Verify spectrum:")
		print(spectrum)

	# Calculate initial exposure without attenuator
	spectrum["Exposure 1"] =\
		spectrum["Emid [keV]"] * \
		spectrum["S(E)"] * \
		spectrum["mu_en/rho (air) (cm^2/g)"] * \
		spectrum["dE [keV]"]
	exposure1 = spectrum["Exposure 1"].sum()

	# get rid of overflow, invalid value warnings
	warnings.filterwarnings("ignore")
	spectrum["Exposure 2"] =\
		spectrum["Emid [keV]"] * \
		spectrum["S(E)"] * \
		spectrum["mu_en/rho (air) (cm^2/g)"] * \
		spectrum["dE [keV]"] * \
		np.exp(
			-spectrum[f"mu ({att_material_name}) (cm-1)"] * thick*1e-1)  # in mm!
	exposure2 = spectrum["Exposure 2"].sum()
	return exposure2 / exposure1


def e_mean(spectrum):
	"""
	A function to calculate the fluence-weighted mean energy of the spectrum

	:param spectrum: DataFrame in a format:
		| Emin [keV] | Emax [keV] | Emid [keV] | dE [keV] | S(E) |
	:type spectrum: pd.DataFrame
	:return: the fluence-weighted mean energy of spectrum
	:rtype: float
	"""
	numerator =\
		(spectrum["Emid [keV]"] * spectrum["S(E)"] * spectrum["dE [keV]"]).sum()
	denominator =\
		(spectrum["S(E)"] * spectrum["dE [keV]"]).sum()

	return numerator/denominator


def e_eff(
		spectrum, att_material="z13", att_material_rho=2.699,
		mu_source="nist", noprint=True):
	"""
	A function to calculate the fluence-weighted mean energy of spectrum

	:param spectrum: DataFrame in a format:
		| Emin [keV] | Emax [keV] | Emid [keV] | dE [keV] | S(E) |
	:type spectrum: pd.DataFrame
	:param att_material: name of the attenuation material from NIST database,
		could be also specified by fraction in dict {Z:frac, ...}
	:type att_material: str | dict
	:param att_material_rho: density of attenuation material
	:type att_material_rho: float
	:param mu_source: "nist" or "penelope"
	:type mu_source: str
	:param noprint: print information while executing or not
	:type noprint: bool
	:return: the fluence-weighted mean energy of spectrum
	:rtype: float
	"""
	def cost_function(energy):
		"""
		A cost function is a measure of how wrong the model is in terms of its
		ability to estimate the relationship between X and y.
		Here MSE(Mean Squared Error) is used as a metric of how close our value to
		the desired.

		:param energy: energy of the photon beam is a free parameter
		:type energy: float
		:return: Mean Squared Error
		:rtype: float
		"""
		mse = np.power(mu - mu_att(energy), 2)
		if not noprint:
			print(mse)
		return mse

	mu = np.log(2.0) / xvl(
		spectrum=spectrum, att_material=att_material,
		att_material_rho=att_material_rho, red_fraction=0.5,
		mu_source=mu_source, noprint=noprint) * 10  # convert to cm
	mu_att = GetMu(
		medium=att_material, rho=att_material_rho,
		mu_type="mass", mu_source=mu_source)

	# TODO: upper limit is hardcoded to the 1000 keV
	return minimize_scalar(
		cost_function, bounds=(1, 1000), method="bounded").x


def fit_exp(x, a, b, c):
	return a*np.exp(-b * x)+c


# TODO: Improve Penelope mu data
class GetMu:
	def __init__(self, medium, rho: float, mu_type="mass", mu_source="nist"):
		"""
		Class returns interp1d object:
		function returns coefficient for specified photon energy (in keV)
		for X-Ray attenuation coefficients (mass or energy) which were loaded
		from NIST or Penelope databases. Only "z13", "z29" and "air" materials
		are available for Penelope case!

		:param medium: required material in dict {Z:frac, ...} or by name
		:type medium: dict, str
		:param rho: density of material in g/cm^3, if rho=1.0, mu/rho will be
			returned
		:type rho: float
		:param mu_type: type of attenuation coefficients "mass" or "energy"
		:type mu_type: str
		:param mu_source: "nist" or "penelope"
		:type mu_source: str
		:return: Interpolator object
		:rtype: scipy.interpolate.interp1d
		"""
		self.medium = medium
		self.rho = rho
		self.mu_type = mu_type
		self.mu_source = mu_source
		names = ["Energy (MeV)", "mu/rho (cm^2/g)", "mu_en/rho (cm^2/g)"]

		def resolve_df(path):
			df = pd.read_csv(
				path,
				sep='\t', skiprows=1,
				names=names, usecols=names, header=None)
			x_dat = df["Energy (MeV)"].to_numpy() * 1e3  # convert to keV
			if mu_type == "energy":
				y_dat = df["mu_en/rho (cm^2/g)"].to_numpy()*self.rho  # convert to mu
			else:
				y_dat = df["mu/rho (cm^2/g)"].to_numpy()*self.rho  # convert to mu
			return np.log(x_dat), np.log(y_dat)  # !!! LOG SCALE !!!

		if self.mu_source == "nist":
			if isinstance(self.medium, dict):
				self.interpolators = {}
				for key in self.medium.keys():
					self.interpolators[key] = \
						interpolate.interp1d(
							*resolve_df(f"data/nist/z{key:02d}.tsv"))
			else:
				self.interpolator = interpolate.interp1d(
					*resolve_df(f"data/nist/{medium}.tsv"))
		elif self.mu_source == "pene":
			if self.mu_type == "energy":  # TODO: only air for now!
				x = \
					np.array(
						load_json(
							f"data/pene/pene_muen_air.dat")["photon energy"])*1e3
				y = \
					np.array(load_json(
						f"data/pene/pene_muen_air.dat")["muen_over_rho_air"])*self.rho
				self.interpolator = interpolate.interp1d(np.log(x), np.log(y))
			else:  # TODO: get rid of this dict hack
				__ = {
					key: value for key, value in zip(
						[f"z{i:02d}" for i in range(1, 93)], [i-1 for i in range(1, 93)])}
				x = \
					np.array(
						load_json(
							f"data/pene/pene_mu.dat")["photon energy"][
							__[self.medium]])*1e3
				y = \
					np.array(load_json(
						f"data/pene/pene_mu.dat")["mu_over_rho"][
						__[self.medium]])*self.rho
				self.interpolator = interpolate.interp1d(np.log(x), np.log(y))
		else:
			exit("No mu data!")

	def __call__(self, x):
		if isinstance(self.medium, dict):
			result = 0
			for key in self.interpolators.keys():
				result +=\
					np.exp(self.interpolators[key](np.log(x))) * \
					self.medium[key]
		else:
			result = np.exp(self.interpolator(np.log(x)))  # !!! LOG !!!
		return result


# TODO: plot graph thick vs. transmission
# TODO: make LaTeX report method to save results into pdf
class SpekQual:
	def __init__(
			self, filename, data_format="fluka", fluka_det_n=0,
			fluka_row_drop=2, fluka_save=False, spekpy_sep="\t",
			spekpy_char=False, mu_source="nist", noprint=True):
		"""
		Class returns spectrum quality object with all the characteristics
		for spectrum

		:param filename: name of the file to load spectrum from
		:type filename: str
		:param data_format: format of the spectrum file "fluka", "spekpy"
		:type data_format: str
		:param fluka_det_n: only for the FLUKA data format option:
			detector number in file to load (Detector n)
		:type fluka_det_n: int
		:param fluka_row_drop: only for the FLUKA spectrum data format option:
			how many rows to drop from the beginning (cut to 1 keV)
		:type fluka_row_drop: int
		:param fluka_save: only for the FLUKA spectrum data format option:
			save extracted spectrum data from FLUKA file to the "result"
		:type fluka_save: bool
		:param spekpy_sep: only for the SpekPy spectrum data format option:
			delimiter used to separate columns in file
		:type spekpy_sep: str
		:param spekpy_char: only for the SpekPy spectrum data format option:
			flag to get column with characteristic fluence otherwise
			total fluence will be considered
		:type spekpy_char: bool
		:param mu_source: "nist" or "penelope"
		:type mu_source: str
		:param noprint: print debug information while executing or not
		:type noprint: bool

		:return: object representing beam spectrum qualities
		"""
		self.filename = filename
		if data_format == "fluka":
			self.spectrum = load_fluka_data(
				filename, det_n=fluka_det_n, row_drop=fluka_row_drop,
				save=fluka_save, noprint=noprint)
		elif data_format == "spekpy":
			self.spectrum = load_spekpy_data(
				filename, sep=spekpy_sep, char=spekpy_char,
				noprint=noprint)
		else:
			exit("Unknown format of spectrum data!")

		self.mu_source = mu_source

		self.hvl1_al, self.hvl2_al, self.qvl_al, self.tvl_al = xvl(
			self.spectrum, att_material="z13", att_material_rho=2.699,
			mu_source=self.mu_source, noprint=noprint)

		self.hvl1_cu, self.hvl2_cu, self.qvl_cu, self.tvl_cu = xvl(
			self.spectrum, att_material="z29", att_material_rho=8.96,
			mu_source=self.mu_source, noprint=noprint)

		# Calculate the homogeneity index hi , which provides a sense
		# of the spectral width and is unity for mono-energetic photons.
		self.hi_al = self.hvl1_al/self.hvl2_al
		self.hi_cu = self.hvl1_cu/self.hvl2_cu

		# The fluence-weighted mean energy of spectrum
		self.emean = e_mean(spectrum=self.spectrum)
		self.eeff_al = e_eff(
			spectrum=self.spectrum, mu_source=self.mu_source, noprint=noprint)
		self.eeff_cu = e_eff(
			self.spectrum, att_material="z29",
			att_material_rho=8.96, mu_source=self.mu_source, noprint=noprint)

		self.__results =\
			f"spek_qual total beam quality characteristics " \
			f"for the spectrum {self.filename}\n\n" \
			f"HVL1 (Al / Cu):\t{self.hvl1_al:.6} / {self.hvl1_cu:.6} mm\n" \
			f"HVL2 (Al / Cu):\t{self.hvl2_al:.6} / {self.hvl2_cu:.6} mm\n" \
			f"QVL (Al / Cu):\t{self.qvl_al:.6} / {self.qvl_cu:.6} mm\n" \
			f"TVL (Al / Cu):\t{self.tvl_al:.6} / {self.tvl_cu:.6} mm\n" \
			f"hi (Al / Cu):\t{self.hi_al:.6} / {self.hi_cu:.6} mm\n" \
			f"Eeff (Al / Cu):\t{self.eeff_al:.6} / {self.eeff_cu:.6} keV\n" \
			f"Emean:\t{self.emean:.6} keV"

	def layer_al(self, fraction):
		"""
		Calculate layer of Al to reduce air Kerma value to the specified
		fraction value

		:param fraction: fraction to reduce to
		:type fraction: float
		:return: layer thickness of Al in mm
		:rtype: float
		"""
		return xvl(
			spectrum=self.spectrum, att_material="z13", att_material_rho=2.699,
			red_fraction=fraction, mu_source=self.mu_source)

	def layer_cu(self, fraction):
		"""
		Calculate layer of Cu to reduce air Kerma value to the specified
		fraction value

		:param fraction: fraction to reduce to
		:type fraction: float
		:return: layer thickness of Al in mm
		:rtype: float
		"""
		return xvl(
			spectrum=self.spectrum, att_material="z29", att_material_rho=8.96,
			red_fraction=fraction, mu_source=self.mu_source)

	def layer(self, fraction, att_material="z13", att_material_rho=2.699):
		"""
		Calculate layer of the specified material to reduce air Kerma value to
		the specified fraction value

		:param fraction: fraction to reduce to
		:type fraction: float
		:param att_material: material to use
		:type att_material: str
		:param att_material_rho: density of material
		:type att_material_rho: float
		:return: layer thickness of material
		:rtype: float
		"""
		return xvl(
			spectrum=self.spectrum, att_material=att_material,
			att_material_rho=att_material_rho, red_fraction=fraction,
			mu_source=self.mu_source)

	def print_all(self):
		"""
		Print all the quality characteristics for the spectrum
		"""
		print(self.__results)

	def save_all(self):
		"""
		Save all the quality characteristics into file using filename
		"""
		with open(
				Path.joinpath(
					Path("result"),
					Path(self.filename).stem+"_spekqual.dat"), mode="w") as f:
			f.write(
				f"spek_qual beam quality characteristics "
				f"for the spectrum {self.filename}\n\n"
				f"HVL1 (Al / Cu):\t{self.hvl1_al:.6} / {self.hvl1_cu:.6} mm\n"
				f"HVL2 (Al / Cu):\t{self.hvl2_al:.6} / {self.hvl2_cu:.6} mm\n"
				f"QVL (Al / Cu):\t{self.qvl_al:.6} / {self.qvl_cu:.6} mm\n"
				f"TVL (Al / Cu):\t{self.tvl_al:.6} / {self.tvl_cu:.6} mm\n"
				f"hi (Al / Cu):\t{self.hi_al:.6} / {self.hi_cu:.6} mm\n"
				f"Eeff (Al / Cu):\t{self.eeff_al:.6} / {self.eeff_cu:.6} keV\n"
				f"Emean:\t{self.emean:.6} keV")

	def get_trans_curve(
			self, att_material="z13", att_material_rho=2.699, step=1e-2):
		"""
		Attenuate spectrum on a specified material thickness and return
		attenuation curve (x, y values)

		:param att_material: name of the attenuation material from NIST database,
			could be also specified by fraction in dict {Z:frac, ...}
		:type att_material: str | dict
		:param att_material_rho: density of attenuation material
		:type att_material_rho: float
		:param step: step for teh thickness calculation
		:type step: float
		:return: result as fraction on specified thickness
		:rtype: float
		"""
		x, y = np.array([]), np.array([])

		thick = 0.0
		frac = 1.0

		while frac > 0.1:  # attenuate till 10%
			x = np.append(x, thick)
			frac = attenuate(
				spectrum=self.spectrum, thick=thick,
				att_material=att_material, att_material_rho=att_material_rho,
				mu_source=self.mu_source)
			y = np.append(y, frac)
			thick += step
		return x, y


def get_args():
	"""
	Get command-line arguments
	"""
	descrp =\
		'spek_qual script calculates beam quality characteristics for the' \
		' specified X-ray spectrum'
	parser = argparse.ArgumentParser(
		description=descrp,
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument(
		'input',
		metavar='spec.dat',
		help='Input spectrum file')

	parser.add_argument(
		'-f',
		'--format',
		help='Format of the spectrum file data: "fluka", "spekpy"',
		metavar='format',
		type=str,
		default='fluka')

	parser.add_argument(
		'-fdn',
		'--fluka_det_n',
		help='Only for the  spectrum data format: detector number in file to '
		'load (Detector n)',
		metavar='fluka_det_n',
		type=int,
		default=1)

	parser.add_argument(
		'-frd',
		'--fluka_row_drop',
		help='Only for the FLUKA spectrum data format: how many rows to drop from '
		'beginning of the file (cut to 1 keV)',
		metavar='fluka_row_drop',
		type=int,
		default=0)

	parser.add_argument(
		'-fs',
		'--fluka_save',
		help='Only for the FLUKA spectrum data format: '
		'save the extracted spectrum data from FLUKA file to the "result" foolder',
		action='store_true')

	parser.add_argument(
		'-ss',
		'--spekpy_sep',
		help='Only for the SpekPy spectrum data format: '
		'delimiter used to separate columns in file with spectrum',
		metavar='spekpy_sep',
		type=str,
		default="\t")

	parser.add_argument(
		'-ms',
		'--mu_source',
		help='Choose mu data: "nist" or "pene"',
		metavar='mu_source',
		type=str,
		default="nist")

	parser.add_argument(
		'-sc',
		'--spekpy_char',
		help='Only for the SpekPy spectrum data format: '
		'flag to get column with the characteristic fluence, otherwise'
		' total fluence will be used for the beam qualification',
		action='store_true')

	parser.add_argument(
		'-v',
		'--verbose',
		help='Print information while execution',
		action='store_true')

	args = parser.parse_args()

	if args.mu_source != "nist" and args.mu_source != "pene":
		parser.error('Only "nist" or "pene" mu data are currently available!')

	if args.format != "fluka" and args.format != "spekpy":
		parser.error('Only "fluka" or "spekpy" formats are currently supported!')

	# Check if input file exists, rise error if not
	if os.path.isfile(f"data/spectra/{args.input}"):
		if args.verbose:
			print(f'Input spectrum file "{args.input}" exists!')
	else:
		parser.error(
			f'Input spectrum file "{args.input}" does not exist!\n'
			'Please put this file into "./data/spectra" directory!')

	# Check numbers
	if args.fluka_det_n < 0:
		parser.error(f'Number of detector must be greater than 0!')
	if args.fluka_row_drop < 0:
		parser.error(f'Number of rows to drop must be non negative!')

	return args


def test():
	return None


def main():

	args = get_args()  # receive args

	# Get SpekQual object
	quality = SpekQual(
		filename=args.input, data_format=args.format,
		fluka_det_n=args.fluka_det_n,
		fluka_row_drop=args.fluka_row_drop,
		fluka_save=args.fluka_save, spekpy_sep=args.spekpy_sep,
		spekpy_char=args.spekpy_char, mu_source=args.mu_source,
		noprint=not args.verbose)

	quality.print_all()
	quality.save_all()


if __name__ == '__main__':
	main()
