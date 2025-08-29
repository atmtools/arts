{% if objname != 'pyarts3.workspace.Workspace' %}
{% if not module.startswith('pyarts3.arts') %}

{{ name }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   .. automethod:: __init__
      :noindex:

   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree:
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
      :toctree:
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% endif %}
{% endif %}
