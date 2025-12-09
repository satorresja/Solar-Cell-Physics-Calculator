# Solar Cell Physics Calculator

A Python-based calculator for modeling heterojunction solar cell performance. This tool computes current-voltage (I-V) characteristics, efficiency, fill factor, and other key parameters based on material properties and spectral data.

## Overview

This software simulates the behavior of a two-layer heterojunction solar cell by calculating:
- Current density contributions from different regions (emitter, base, space charge region)
- Dark current based on recombination mechanisms
- Cell performance metrics (efficiency, fill factor, Voc, Jsc)

## Requirements

- Python 3.7 or higher
- Required packages:
  ```bash
  pip install numpy pandas scipy
  ```

## Quick Start

1. **Prepare your input files** (CSV format):
   - `REFLECTANCE.csv` - Spectral reflectance data
   - `Transmitancia.csv` - Spectral transmittance data
   - `ESAM.csv` - Solar spectrum photon flux (AM1.5 standard)

2. **Run the calculator**:
   ```bash
   python solar_cell.py
   ```

3. **View results** in `solar_cell_results.txt`

## Example Data

This repository includes example datasets for testing:
- **REFLECTANCE.csv**: Sample reflectance spectrum
- **Transmitancia.csv**: Sample transmittance spectrum
- **ESAM.csv**: AM1.5 solar spectrum data

The AM1.5 solar spectrum is the standard reference defined by NREL:  
https://www.nrel.gov/grid/solar-resource/spectra-am1.5

## Input File Format

All input files should be CSV format with two columns (no headers):

### Reflectance and Transmittance Files
```
wavelength_1, value_1
wavelength_2, value_2
...
```

**Example** (`REFLECTANCE.csv`):
```
280, 0.057784785
280.5, 0.057529785
281, 0.057276816
...
```

### Solar Spectrum File (ESAM)
```
wavelength_1, photon_flux_1
wavelength_2, photon_flux_2
...
```

**Note**: The software uses the **second column** as photon flux values.

## Customizing Material Properties

You can modify material properties by editing the `MaterialProperties` class in `solar_cell.py`:

### Key Parameters You Can Change

#### **Band Gaps** (in eV)
```python
band_gap_p: float = 1.17  # p-layer band gap
band_gap_n: float = 2.42  # n-layer band gap
```
Change these to match your materials (e.g., Si = 1.12 eV, GaAs = 1.42 eV)

#### **Doping Concentrations** (in cm⁻³)
```python
doping_acceptor: float = 2e16  # Acceptor doping (p-layer)
doping_donor: float = 1e17     # Donor doping (n-layer)
```
Typical range: 10¹⁴ to 10¹⁹ cm⁻³

#### **Layer Thicknesses** (in cm)
```python
width_n: float = 0.1e-4   # n-layer thickness (1 μm)
width_p: float = 2e-4     # p-layer thickness (20 μm)
```

#### **Diffusion Parameters**
```python
diffusion_length_p: float = 2.9e-6      # Minority carrier diffusion length in n-layer
diffusion_length_n: float = 2.31e-4     # Minority carrier diffusion length in p-layer
diffusion_coefficient_p: float = 0.65   # Diffusion coefficient (cm²/s)
diffusion_coefficient_n: float = 1.05   # Diffusion coefficient (cm²/s)
```

#### **Surface Recombination Velocities** (in cm/s)
```python
surface_recombination_p: float = 1e7  # Front surface
surface_recombination_n: float = 1e1  # Back surface
```

#### **Other Parameters**
```python
temperature: float = 300.0              # Temperature (K)
epsilon_p: float = 13.6                 # Relative permittivity (p-layer)
epsilon_n: float = 10.0                 # Relative permittivity (n-layer)
electron_affinity_p: float = 4.235      # Electron affinity (eV)
electron_affinity_n: float = 4.3        # Electron affinity (eV)
absorption_edge_index: int = 900        # Wavelength index for absorption limit
```

## How to Modify Parameters

### Method 1: Edit the Default Values
Open `solar_cell.py` and locate the `MaterialProperties` class (around line 18). Change the default values directly:

```python
@dataclass
class MaterialProperties:
    temperature: float = 300.0
    band_gap_p: float = 1.42    # Changed from 1.17 to 1.42 for GaAs
    band_gap_n: float = 3.4     # Changed from 2.42 to 3.4 for ZnO
    # ... modify other parameters as needed
```

### Method 2: Create Custom Properties in Code
If you're comfortable with Python, create a custom configuration:

```python
from solar_cell import MaterialProperties, SolarCellCalculator

# Create custom material properties
my_materials = MaterialProperties(
    band_gap_p=1.42,        # GaAs band gap
    band_gap_n=3.4,         # ZnO band gap
    doping_acceptor=5e16,   # Custom doping
    doping_donor=2e17       # Custom doping
)

# Run calculator with custom properties
calculator = SolarCellCalculator(my_materials)
voltage, current = calculator.calculate_iv_characteristics(
    'REFLECTANCE.csv',
    'Transmitancia.csv',
    'ESAM.csv'
)
```

## Output

The software generates `solar_cell_results.txt` with:
```
Efficiency: XX.XX
Fill Factor: XX.XX%
Short Circuit Current: XX.XX mA/cm²
Open Circuit Voltage: X.XXXX V
Maximum Power: XX.XX
```

### Output Metrics Explained

- **Efficiency**: Solar cell power conversion efficiency (in the units of max power)
- **Fill Factor (FF)**: Ratio of maximum power to the product of Voc and Jsc (%)
- **Short Circuit Current (Jsc)**: Current when voltage is zero (mA/cm²)
- **Open Circuit Voltage (Voc)**: Voltage when current is zero (V)
- **Maximum Power**: Maximum power output point

## Troubleshooting

### Common Issues

**Problem**: "File not found" error  
**Solution**: Ensure your CSV files are in the same directory as `solar_cell.py`

**Problem**: "No module named 'numpy'" error  
**Solution**: Install required packages: `pip install numpy pandas scipy`

**Problem**: Unrealistic results (negative values, very high/low efficiency)  
**Solution**: Check your material parameters. Common issues:
- Band gaps should be in eV (typical range: 0.5-4.0)
- Doping should be in cm⁻³ (typical range: 10¹⁴-10¹⁹)
- Layer thicknesses in cm (typical range: 10⁻⁶ to 10⁻² cm)

**Problem**: CSV parsing errors  
**Solution**: Ensure your CSV files:
- Use commas as separators
- Have no headers
- Have two columns: wavelength and value
- Use decimal points (not commas) for numbers

## Understanding the Physics

The calculator implements a comprehensive model including:

1. **Generation current** from photon absorption in:
   - Emitter (n-layer)
   - Base (p-layer)
   - Space charge region (depletion zone)
   - Contributions from reflected light

2. **Dark current** from:
   - Diffusion of minority carriers
   - Recombination in the space charge region

3. **Key physics equations**:
   - Built-in potential from band alignment
   - Depletion width calculation
   - Absorption coefficient (square root law)
   - Minority carrier transport equations

## Tips for Best Results

1. **Wavelength range**: Ensure your spectral data covers 280-1300 nm
2. **Data resolution**: Use 0.5-1 nm wavelength steps for accuracy
3. **Consistent units**: All three CSV files should use the same wavelength scale
4. **Material matching**: Ensure band gaps, electron affinities, and permittivities match your actual materials
5. **Physical constraints**: Verify that diffusion lengths are smaller than layer thicknesses.

## Citation

If you use this software in your research, please cite the underlying physical models and the AM1.5 solar spectrum:

- AM1.5 Reference Spectrum: ASTM G173-03
- NREL Solar Spectra: https://www.nrel.gov/grid/solar-resource/spectra-am1.5
- https://www.elsevier.es/en-revista-journal-applied-research-technology-jart-81-articulo-design-thin-film-solar-cells-S1665642317300950
- https://www.sciencedirect.com/science/article/abs/pii/S0030402618312749?via%3Dihub

## License

This software is provided as-is for educational and research purposes.

## Support

For questions or issues:
1. Check that all CSV files are properly formatted
2. Verify material parameters are physically reasonable
3. Ensure Python packages are correctly installed

## Version History

- **v1.0**: Initial release with heterojunction solar cell modeling
  - Added example datasets (transmittance, reflectance, AM1.5 spectrum)
  - Based on NREL AM1.5 reference spectrum
