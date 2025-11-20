# HiM
​​HiM (Hybrid iterative ptychographic super-resolution Microscopy)​​ is an advanced computational imaging method that synergistically combines a honeycomb array detection architecture with a hybrid high-frequency expansion algorithm. It achieves ~60nm super-resolution while significantly reducing photobleaching effects. This repository contains implementation code for super-resolution reconstruction and sample test data.

## Quick Start
1. Download and extract `data.7z`
2. Open `HiM_main.m` in MATLAB
3. Check and modify `flagGpu` parameter if needed:
   - Set `flagGpu = 0` for CPU-only execution
   - Set `flagGpu = 1` for GPU-accelerated execution (default)
4. Execute `HiM_main.m` to perform the reconstruction

## Enviroment
- **Operating System**: Windows 10
- **Software**: MATLAB R2020a or later
- **Hardware Requirements**:
  - **Minimum**: CPU only
  - **Recommended**: NVIDIA GPU with CUDA support
  - **Tested Configuration**: NVIDIA Quadro P4000
