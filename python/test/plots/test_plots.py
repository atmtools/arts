# %%
from os.path import join, dirname
import matplotlib.pyplot as plt
import pyarts3.xml as axml
from pyarts3.plots import plot_arts_lookup
from pyarts3.arts import ArrayOfArrayOfSpeciesTag, SpeciesTag


class TestPlots:
    """Testing the ARTS plotting functions."""
    def test_plot_lookup(self):
        lookup_file = join(dirname(__file__), 'reference',
                           'abs_lookup_small.xml')
        fig, ax = plot_arts_lookup(axml.load(lookup_file))
        
        fig.suptitle('Lookup table opacities')
        fig.subplots_adjust(top=0.88)
        plt.savefig("lookup_opacities.pdf")
        
        lookup_file = join(dirname(__file__), 'reference',
                           'abs_lookup_small.xml')
        fig, ax = plot_arts_lookup(
            axml.load(lookup_file),
            species=ArrayOfArrayOfSpeciesTag([[SpeciesTag("N2O")],
                                              [SpeciesTag("O3")]]),
            opacity=False)
        
        fig.suptitle('Lookup table absorption cross sections [m$^2$]')
        fig.subplots_adjust(top=0.88)
        plt.savefig("lookup_crosssections.pdf")

if __name__ == "__main__":
    x = TestPlots()
    x.test_plot_lookup()
