function [Ex,Ey,Ez,varargout] = singleobjectivepsf(Obj,Beam,Scope,Option)
%SINGLEOBJECTIVEPSF calculates the PSF of single objective
%
% *************************************************************************
% originates from Dr. Hao
% optimized by LIU Xin
% liuxin2018@zju.edu.cn
% Apr.24, 2021

if ~isfield(Option,'FocusingMethod')
    Option.FocusingMethod = 'mm';
    warning(['FocusingMethod is set as default: ', Option.FocusingMethod]);
end

if ~isfield(Option,'UseGpu')
    Option.UseGpu = 0;
    warning(['UseGpu is set as default: ', num2str(Option.UseGpu)]);
end

if ~isfield(Option,'Precision')
    Option.Precision = 'double';
    warning(['Precision is set as default: ', Option.Precision]);
end

if ~isfield(Option,'BatchSize')
    Option.BatchSize = 1;
    warning(['BatchSize is set as default: ', num2str(Option.BatchSize)]);
end

if ~isfield(Option,'TimeIt')
    Option.TimeIt = 0;
    warning(['TimeIt is set as default: ', num2str(Option.TimeIt)]);
end

switch Option.FocusingMethod
    case 'int'  % direct integral (most accurate and flexible, but slowest)
        [Ex,Ey,Ez] = singleobjectivepsf_noT_INT(Obj,Beam,Scope,Option);
    case 'mm'  % matrix multiplication (most accurate and fastest, middling flexible)
        [Ex,Ey,Ez] = singleobjectivepsf_noT_MM(Obj,Beam,Scope,Option);
    case 'czt'  % 2D czt (most accurate and fastest, middling flexible)
        [Ex,Ey,Ez] = singleobjectivepsf_noT_CZT(Obj,Beam,Scope,Option);
    case 'czt3'  % 3D czt (middling accurate and fast, middling flexible)
        [Ex,Ey,Ez] = singleobjectivepsf_noT_CZT_3D(Obj,Beam,Scope,Option);
    case 'fft'  % 2D fft (most accurate, but very slow due to zeros padding, least flexible)
        [Ex,Ey,Ez,scopeNew] = singleobjectivepsf_noT_FFT(Obj,Beam,Scope,Option);
        varargout{1} = scopeNew;
    case 'nufft'  % 2D nufft (most accurate, middling flexible and fast)
        [Ex,Ey,Ez] = singleobjectivepsf_noT_NUFFT(Obj,Beam,Scope,Option);
    otherwise
        error('No such a mode, please try again!');
end

[Ex,Ey,Ez] = tranelectricfield(Scope,Ex,Ey,Ez);  % (unit: V/m)
end