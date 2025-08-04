function session_name_str = para_to_namestr(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior)       
b_PF_str                = sprintf('%.2f', b_PF);
b_PF_str                = strrep(b_PF_str, '.', '_');

cardinal_delta_str      = sprintf('%.2f',cardinal_delta);
cardinal_delta_str      = strrep(cardinal_delta_str, '.', '_');

cardinal_prior_str      = sprintf('%.2f',cardinal_prior);
cardinal_prior_str      = strrep(cardinal_prior_str, '.', '_');

oblique_delta_str       = sprintf('%.2f',oblique_delta);
oblique_delta_str       = strrep(oblique_delta_str, '.', '_');

oblique_prior_str       = sprintf('%.2f',oblique_prior);
oblique_prior_str       = strrep(oblique_prior_str, '.', '_');


session_name_str = sprintf('bPF_%s_cardinal_delta_%s_prior_%s_oblique_delta_%s_prior_%s.mat',...
                    b_PF_str, cardinal_delta_str, cardinal_prior_str, oblique_delta_str, oblique_prior_str);
end