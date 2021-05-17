function [Ve_freq_octave, Ve_STL_octave] = fun_octave(Ma)

Winc = 9.5198e-4;
Ve_Wtrans_narrow = 10.^(log10(Winc) - Ma(:,2)./10);
[Ve_freq_octave, Ve_Wtrans_octave] = fun_narrow_to_one_third_octave(Ma(:,1),Ve_Wtrans_narrow);
Ve_STL_octave =  10*(log10(Winc) - log10(Ve_Wtrans_octave));

end