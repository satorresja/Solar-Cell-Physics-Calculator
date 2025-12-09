import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Tuple
from scipy.interpolate import interp1d


@dataclass
class PhysicalConstants:
    q = 1.6e-19
    k_boltzmann = 8.617332e-5
    planck_ev = 4.1356e-15
    speed_of_light = 3e17
    epsilon_0 = 8.8541e-14


@dataclass
class MaterialProperties:
    temperature: float = 300.0
    band_gap_p: float = 1.17
    band_gap_n: float = 2.42
    
    diffusion_length_p: float = 2.9e-6
    diffusion_length_n: float = 2.31e-4
    
    surface_recombination_p: float = 1e7
    surface_recombination_n: float = 1e1
    
    diffusion_coefficient_p: float = 0.65
    diffusion_coefficient_n: float = 1.05
    
    width_n: float = 0.1e-4
    width_p: float = 2e-4
    
    epsilon_p: float = 13.6
    epsilon_n: float = 10.0
    
    doping_acceptor: float = 2e16
    doping_donor: float = 1e17
    
    density_states_conduction_p: float = 2.2e18
    density_states_valence_n: float = 1.8e19
    density_states_conduction_n: float = 2.2e18
    density_states_valence_p: float = 1.8e19
    
    electron_affinity_p: float = 4.235
    electron_affinity_n: float = 4.3
    
    absorption_edge_index: int = 900


class SolarCellCalculator:
    def __init__(self, material_props: MaterialProperties):
        self.props = material_props
        self.const = PhysicalConstants()
        
    def load_spectral_data(self, reflectance_file: str, transmittance_file: str, 
                          photon_flux_file: str = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        reflectance_data = pd.read_csv(reflectance_file, header=None, names=['wavelength', 'reflectance'])
        transmittance_data = pd.read_csv(transmittance_file, header=None, names=['wavelength', 'transmittance'])
        
        wavelength = reflectance_data['wavelength'].values
        reflectance = reflectance_data['reflectance'].values
        transmittance = transmittance_data['transmittance'].values
        
        if photon_flux_file:
            photon_flux_data = pd.read_csv(photon_flux_file, header=None, usecols=[0, 1], 
                                          names=['wavelength', 'flux'], skipinitialspace=True)
            photon_flux_data = photon_flux_data.dropna()
            photon_flux = photon_flux_data['flux'].values
        else:
            photon_flux = np.ones_like(wavelength)
            
        return wavelength, reflectance, transmittance, photon_flux
    
    def calculate_photon_energy(self, wavelength: np.ndarray) -> np.ndarray:
        return self.const.planck_ev * self.const.speed_of_light / wavelength
    
    def calculate_intrinsic_carrier_concentrations(self) -> Tuple[float, float]:
        ni_p = np.sqrt(self.props.density_states_conduction_p * self.props.density_states_valence_p) * \
               np.exp(-self.props.band_gap_p / (2 * self.const.k_boltzmann * self.props.temperature))
        
        ni_n = np.sqrt(self.props.density_states_conduction_n * self.props.density_states_valence_n) * \
               np.exp(-self.props.band_gap_n / (2 * self.const.k_boltzmann * self.props.temperature))
        
        return ni_p, ni_n
    
    def calculate_built_in_potential(self, ni_p: float, ni_n: float) -> float:
        term1 = self.const.k_boltzmann * self.props.temperature * \
                np.log((self.props.doping_acceptor * self.props.doping_donor) / (ni_p * ni_n))
        
        term2 = self.props.electron_affinity_p - self.props.electron_affinity_n
        
        term3 = (self.props.band_gap_p - self.props.band_gap_n) / 2
        
        term4 = -(self.const.k_boltzmann * self.props.temperature / 2) * \
                np.log((self.props.density_states_valence_p * self.props.density_states_conduction_n) / 
                       (self.props.density_states_conduction_p * self.props.density_states_valence_n))
        
        return term1 + term2 + term3 + term4
    
    def calculate_depletion_widths(self, voltage: np.ndarray, built_in_potential: float) -> Tuple[np.ndarray, np.ndarray]:
        factor = 2 * self.props.epsilon_p * self.props.epsilon_n * self.const.epsilon_0
        denominator = self.const.q * (self.props.epsilon_n * self.props.doping_donor + 
                                      self.props.epsilon_p * self.props.doping_acceptor)
        
        depletion_n = np.sqrt((factor * self.props.doping_acceptor * (built_in_potential - voltage)) / 
                             (self.props.doping_donor * denominator))
        
        depletion_p = np.sqrt((factor * self.props.doping_donor * (built_in_potential - voltage)) / 
                             (self.props.doping_acceptor * denominator))
        
        return depletion_n, depletion_p
    
    def calculate_absorption_coefficients(self, photon_energy: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        absorption_constant = 1e5
        
        alpha_n = np.where(photon_energy > self.props.band_gap_n,
                          absorption_constant * np.sqrt(photon_energy - self.props.band_gap_n),
                          0.0)
        
        alpha_p = np.where(photon_energy > self.props.band_gap_p,
                          absorption_constant * np.sqrt(photon_energy - self.props.band_gap_p),
                          0.0)
        
        return alpha_n, alpha_p
    
    def calculate_current_density_emitter(self, photon_flux: np.ndarray, reflectance: np.ndarray,
                                         transmittance: np.ndarray, alpha_n: np.ndarray,
                                         depletion_n: float) -> np.ndarray:
        L_p = self.props.diffusion_length_p
        S_p = self.props.surface_recombination_p
        D_p = self.props.diffusion_coefficient_p
        W_n = self.props.width_n
        
        prefactor = (self.const.q * photon_flux * (1 - reflectance) * transmittance * alpha_n * L_p) / \
                   ((alpha_n * L_p)**2 - 1)
        
        cosh_term = np.cosh((W_n - depletion_n) / L_p)
        sinh_term = np.sinh((W_n - depletion_n) / L_p)
        surface_term = S_p * L_p / D_p
        exp_term = np.exp(-alpha_n * (W_n - depletion_n))
        
        numerator = surface_term + alpha_n * L_p - exp_term * (surface_term * cosh_term + sinh_term)
        denominator = surface_term * sinh_term + cosh_term
        
        return prefactor * (numerator / denominator - alpha_n * L_p * exp_term)
    
    def calculate_current_density_base(self, photon_flux: np.ndarray, reflectance: np.ndarray,
                                      transmittance: np.ndarray, alpha_n: np.ndarray, alpha_p: np.ndarray,
                                      depletion_n: float, depletion_p: float) -> np.ndarray:
        L_n = self.props.diffusion_length_n
        S_n = self.props.surface_recombination_n
        D_n = self.props.diffusion_coefficient_n
        W_n = self.props.width_n
        W_p = self.props.width_p
        
        prefactor = (self.const.q * photon_flux * (1 - reflectance) * transmittance * alpha_p * L_n) / \
                   ((alpha_p * L_n)**2 - 1)
        
        exp_absorption = np.exp(-(alpha_n * W_n + alpha_p * depletion_p))
        surface_term = S_n * L_n / D_n
        
        cosh_term = np.cosh((W_p - depletion_p) / L_n)
        sinh_term = np.sinh((W_p - depletion_p) / L_n)
        exp_width = np.exp(-alpha_p * (W_p - depletion_p))
        
        numerator = alpha_p * L_n - (surface_term * (cosh_term - exp_width) + sinh_term + alpha_p * L_n * exp_width)
        denominator = surface_term * sinh_term + cosh_term
        
        return prefactor * exp_absorption * (numerator / denominator)
    
    def calculate_current_density_space_charge(self, photon_flux: np.ndarray, reflectance: np.ndarray,
                                               transmittance: np.ndarray, alpha_n: np.ndarray, alpha_p: np.ndarray,
                                               depletion_n: float, depletion_p: float) -> np.ndarray:
        exp_n = np.exp(-alpha_n * depletion_n)
        exp_p = np.exp(-alpha_p * depletion_p)
        exp_combined = np.exp(-alpha_n * (self.props.width_n - depletion_n))
        
        return self.const.q * photon_flux * (1 - reflectance) * transmittance * exp_combined * \
               ((1 - exp_n) + exp_n * (1 - exp_p))
    
    def calculate_reflection_contributions(self, photon_flux: np.ndarray, reflectance: np.ndarray,
                                          transmittance: np.ndarray, alpha_n: np.ndarray, alpha_p: np.ndarray,
                                          depletion_n: float, depletion_p: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        W_n = self.props.width_n
        W_p = self.props.width_p
        L_n = self.props.diffusion_length_n
        L_p = self.props.diffusion_length_p
        S_n = self.props.surface_recombination_n
        S_p = self.props.surface_recombination_p
        D_n = self.props.diffusion_coefficient_n
        D_p = self.props.diffusion_coefficient_p
        
        prefactor = self.const.q * photon_flux * (1 - reflectance) * transmittance
        
        exp_abs = np.exp(-(alpha_n * W_n + alpha_p * W_p))
        J_abs = prefactor * (alpha_p * L_n / ((alpha_p * L_n)**2 - 1)) * exp_abs * \
                (alpha_p * L_n - ((S_n * L_n / D_n) * (np.cosh((W_p - depletion_p) / L_n) - 
                 np.exp(-alpha_p * (W_p - depletion_p))) + np.sinh((W_p - depletion_p) / L_n) + 
                 alpha_p * L_n * np.exp(-alpha_p * (W_p - depletion_p))) / 
                ((S_n * L_n / D_n) * np.sinh((W_p - depletion_p) / L_n) + np.cosh((W_p - depletion_p) / L_n)))
        
        exp_vent = np.exp(-(alpha_n * (W_n + depletion_n) + alpha_n * 2 * W_p))
        J_vent = prefactor * (alpha_n * L_p / ((alpha_n * L_p)**2 - 1)) * exp_vent * \
                 (alpha_n * L_p - ((S_p * L_p / D_p) * (np.cosh((W_n - depletion_n) / L_p) - 
                  np.exp(-alpha_n * (W_n - depletion_n))) + np.sinh((W_n - depletion_n) / L_p) + 
                  alpha_n * L_p * np.exp(-alpha_n * (W_n - depletion_n))) / 
                 ((S_p * L_p / D_p) * np.sinh((W_n - depletion_n) / L_p) + np.cosh((W_n - depletion_n) / L_p)))
        
        exp_re_p = np.exp(-(alpha_n * W_n + alpha_p * (2 * W_p - depletion_p)))
        J_re_p = prefactor * exp_re_p * (1 - np.exp(-alpha_p * depletion_p))
        
        exp_re_n = np.exp(-(alpha_n * W_n + alpha_p * 2 * W_p))
        J_re_n = prefactor * exp_re_n * (1 - np.exp(-alpha_n))
        
        return J_abs, J_vent, J_re_p, J_re_n
    
    def calculate_dark_current(self, voltage: np.ndarray, depletion_n: np.ndarray, 
                              depletion_p: np.ndarray, ni_p: float, ni_n: float) -> np.ndarray:
        L_p = self.props.diffusion_length_p
        L_n = self.props.diffusion_length_n
        D_p = self.props.diffusion_coefficient_p
        D_n = self.props.diffusion_coefficient_n
        S_p = self.props.surface_recombination_p
        S_n = self.props.surface_recombination_n
        W_n = self.props.width_n
        W_p = self.props.width_p
        
        tau_p = L_p**2 / D_p
        tau_n = L_n**2 / D_n
        
        n_o = ni_p**2 / self.props.doping_acceptor
        p_o = ni_n**2 / self.props.doping_donor
        
        J_recombination = 1000 * self.const.q * ((depletion_n * ni_n) / tau_p + (depletion_p * ni_p) / tau_n)
        
        cosh_n = np.cosh((W_p - depletion_p) / L_n)
        sinh_n = np.sinh((W_p - depletion_p) / L_n)
        J_o_n = 1000 * (self.const.q * D_n * n_o / L_n) * \
                ((S_n * L_n / D_n) * cosh_n + sinh_n) / ((S_n * L_n / D_n) * sinh_n + cosh_n)
        
        cosh_p = np.cosh((W_n - depletion_n) / L_p)
        sinh_p = np.sinh((W_n - depletion_n) / L_p)
        J_o_p = 1000 * (self.const.q * D_p * p_o / L_p) * \
                ((S_p * L_p / D_p) * cosh_p + sinh_p) / ((S_p * L_p / D_p) * sinh_p + cosh_p)
        
        J_o = J_o_n + J_o_p
        
        thermal_voltage = self.const.k_boltzmann * self.props.temperature
        ideal_factor = np.exp(self.const.q * voltage / thermal_voltage) - 1
        recombination_factor = np.exp(self.const.q * voltage / (2 * thermal_voltage)) - 1
        
        return J_o * ideal_factor + J_recombination * recombination_factor
    
    def calculate_iv_characteristics(self, reflectance_file: str, transmittance_file: str,
                                    photon_flux_file: str = None, num_voltage_points: int = None):
        wavelength, reflectance, transmittance, photon_flux = self.load_spectral_data(
            reflectance_file, transmittance_file, photon_flux_file)
        
        if num_voltage_points is None:
            num_voltage_points = len(wavelength)
        
        photon_energy = self.calculate_photon_energy(wavelength)
        alpha_n, alpha_p = self.calculate_absorption_coefficients(photon_energy)
        ni_p, ni_n = self.calculate_intrinsic_carrier_concentrations()
        built_in_potential = self.calculate_built_in_potential(ni_p, ni_n)
        
        voltage = np.linspace(0, 0.8, num_voltage_points)
        current_density = np.zeros(num_voltage_points)
        
        absorption_limit = self.props.absorption_edge_index
        
        for i, v in enumerate(voltage):
            depletion_n, depletion_p = self.calculate_depletion_widths(np.array([v]), built_in_potential)
            depletion_n = depletion_n[0]
            depletion_p = depletion_p[0]
            
            J_emitter = self.calculate_current_density_emitter(
                photon_flux, reflectance, transmittance, alpha_n, depletion_n)
            
            J_base = self.calculate_current_density_base(
                photon_flux, reflectance, transmittance, alpha_n, alpha_p, depletion_n, depletion_p)
            
            J_scr = self.calculate_current_density_space_charge(
                photon_flux, reflectance, transmittance, alpha_n, alpha_p, depletion_n, depletion_p)
            
            J_abs, J_vent, J_re_p, J_re_n = self.calculate_reflection_contributions(
                photon_flux, reflectance, transmittance, alpha_n, alpha_p, depletion_n, depletion_p)
            
            photocurrent = (np.trapz(J_emitter[:absorption_limit]) + 
                          np.trapz(J_base[:absorption_limit]) +
                          np.trapz(J_scr[:absorption_limit]) +
                          np.trapz(J_abs[:absorption_limit]) +
                          np.trapz(J_vent[:absorption_limit]) +
                          np.trapz(J_re_p[:absorption_limit]) +
                          np.trapz(J_re_n[:absorption_limit]))
            
            dark_current = self.calculate_dark_current(
                np.array([v]), np.array([depletion_n]), np.array([depletion_p]), ni_p, ni_n)[0]
            
            current_density[i] = photocurrent / 10 - dark_current
        
        return voltage, current_density
    
    def calculate_cell_parameters(self, voltage: np.ndarray, current_density: np.ndarray) -> dict:
        positive_mask = current_density > 0
        if np.any(positive_mask):
            zero_crossing_idx = np.where(np.diff(positive_mask.astype(int)) < 0)[0]
            if len(zero_crossing_idx) > 0:
                idx = zero_crossing_idx[0]
                v_oc = interp1d(current_density[idx:idx+2], voltage[idx:idx+2], 
                               fill_value='extrapolate')(0)
            else:
                v_oc = voltage[-1]
        else:
            v_oc = 0.0
        
        j_sc = np.max(current_density)
        
        power = voltage * current_density
        max_power = np.max(power)
        max_power_idx = np.argmax(power)
        v_mp = voltage[max_power_idx]
        j_mp = current_density[max_power_idx]
        
        if v_oc > 0 and j_sc > 0:
            fill_factor = 100 * max_power / (v_oc * j_sc)
        else:
            fill_factor = 0.0
        
        efficiency = max_power
        
        return {
            'efficiency': efficiency,
            'fill_factor': fill_factor,
            'short_circuit_current': j_sc,
            'open_circuit_voltage': v_oc,
            'max_power': max_power,
            'voltage_at_max_power': v_mp,
            'current_at_max_power': j_mp
        }
    
    def save_results(self, results: dict, filename: str = 'solar_cell_results.txt'):
        with open(filename, 'w') as f:
            f.write(f"Efficiency: {results['efficiency']:.2f}\n")
            f.write(f"Fill Factor: {results['fill_factor']:.2f}%\n")
            f.write(f"Short Circuit Current: {results['short_circuit_current']:.2f} mA/cm²\n")
            f.write(f"Open Circuit Voltage: {results['open_circuit_voltage']:.4f} V\n")
            f.write(f"Maximum Power: {results['max_power']:.2f}\n")


if __name__ == "__main__":
    material = MaterialProperties()
    calculator = SolarCellCalculator(material)
    
    voltage, current = calculator.calculate_iv_characteristics(
        'REFLECTANCE.csv',
        'Transmitancia.csv',
        'ESAM.csv'
    )
    
    results = calculator.calculate_cell_parameters(voltage, current)
    calculator.save_results(results)
    
    print("Solar Cell Characteristics:")
    print(f"Efficiency: {results['efficiency']:.2f}")
    print(f"Fill Factor: {results['fill_factor']:.2f}%")
    print(f"Jsc: {results['short_circuit_current']:.2f} mA/cm²")
    print(f"Voc: {results['open_circuit_voltage']:.4f} V")
