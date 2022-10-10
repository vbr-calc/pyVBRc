function VBR = generate_sample_vbr_struct_file(output_dir)
    % generates and saves a VBRc run to use as a lookup table.
    % varies T, phi, dg_um (see genSVranges in this file for full settings) at
    % two frequencies.
    path_to_top_level_vbr = getenv("vbrdir");
    addpath(path_to_top_level_vbr)
    vbr_init

    SVs = genSVranges();
    VBRsettings.ane_meths={'andrade_psp';'xfit_mxw';'eburgers_psp';'xfit_premelt'};
    VBRsettings.freqs=[0.01, 0.1];

    VBR = runVBR(SVs, VBRsettings);
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        save(fullfile(output_dir,'VBRc_sample_LUT.mat'), 'VBR', '-mat-binary')
    else
        save(fullfile(output_dir,'VBRc_sample_LUT.mat'), 'VBR')
    end

end

function [SVs] = genSVranges()
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % [SVs,Ranges] = genSVranges()
  %
  % builds state variable structure and ranges varied.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Ranges.T_K=800:100:1500 + 273;
  Ranges.phi=logspace(-8,-2,10);
  Ranges.dg_um=logspace(-3,-2,5)*1e6;

  Constants.sig_MPa=0.1;
  Constants.P_GPa=2.5;
  Constants.rho=3300;
  Constants.Tsolidus_K=1200+273;

  % get length of each range
  flds=fieldnames(Ranges);
  for ifld=1:numel(flds)
    N.(flds{ifld})=numel(Ranges.(flds{ifld}));
  end

  % build SVs for each var (vectorize this)
  flds=fieldnames(Ranges);
  SVs.T_K = zeros(N.T_K, N.phi, N.dg_um);
  SVs.phi = zeros(N.T_K, N.phi, N.dg_um);
  SVS.dg_um = zeros(N.T_K, N.phi, N.dg_um);
  for iT = 1:N.T_K
      for iphi = 1:N.phi
          for idg = 1:N.dg_um
              SVs.T_K(iT, iphi, idg) = Ranges.T_K(iT);
              SVs.phi(iT, iphi, idg) = Ranges.phi(iphi);
              SVs.dg_um(iT, iphi, idg) = Ranges.dg_um(idg);
          end
      end
  end

  % fill in the other constants
  flds=fieldnames(Constants);
  onz=ones(size(SVs.T_K));
  for ifld=1:numel(flds)
    SVs.(flds{ifld})=Constants.(flds{ifld}) * onz;
  end

end



function VBR = runVBR(SVs,VBRsettings);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Box = runVBR(Box,VBRsettings)
  %
  % generates or loads the Box
  %
  % Parameters
  % ----------
  % SVs          the state variables
  % VBRsettings  structure of some VBR settings
  %
  % Output
  % ------
  % VBR          the VBR structure
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Load and set shared VBR parameters
  VBR.in.elastic.methods_list={'anharmonic','anh_poro'};
  VBR.in.viscous.methods_list={'HK2003'};
  VBR.in.anelastic.methods_list=VBRsettings.ane_meths;
  VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
  VBR.in.elastic.anharmonic.Gu_0_ol = 75.5; % olivine reference shear modulus [GPa]

  VBR.in.SV=SVs;
  VBR.in.SV.f=VBRsettings.freqs;
  disp('Calculating material properties....')
  [VBR] = VBR_spine(VBR) ;
end
