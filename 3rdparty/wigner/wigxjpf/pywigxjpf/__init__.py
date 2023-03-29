try:
    from pywigxjpf import \
        wig_table_init, wig_temp_init, wig_temp_free, wig_table_free, \
        wig3jj, wig6jj, wig9jj, wig3jj_array, wig6jj_array, wig9jj_array
except ImportError:
    from pywigxjpf.pywigxjpf import \
        wig_table_init, wig_temp_init, wig_temp_free, wig_table_free, \
        wig3jj, wig6jj, wig9jj, wig3jj_array, wig6jj_array, wig9jj_array
