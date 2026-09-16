"""Show effective workspace defaults in Sphinx signatures."""

from docutils import nodes
from sphinx import addnodes


def format_workspace_signatures(app, doctree):
    from pyarts3.arts import globals as arts_globals

    methods = arts_globals.workspace_methods()
    for signature in doctree.findall(addnodes.desc_signature):
        if (signature.get('module'), signature.get('class')) not in {
            ('pyarts3.workspace', 'Workspace'),
            ('pyarts3.arts', 'CxxWorkspace'),
        }:
            continue
        method_name = signature.get('fullname', '').rsplit('.', 1)[-1]
        if method_name not in methods:
            continue
        method = methods[method_name]
        workspace_arguments = set(method.input) | set(method.output)
        generic_defaults = dict(zip(method.gin, method.gin_value))
        for parameter in signature.findall(addnodes.desc_parameter):
            name = parameter.children[0].astext()
            if name not in workspace_arguments and (
                    name not in generic_defaults or name in method.gout):
                continue
            defaults = [n for n in parameter.findall(nodes.inline)
                        if 'default_value' in n.get('classes', [])]
            if len(defaults) != 1 or defaults[0].astext() != 'None':
                continue
            default = defaults[0]
            if name in workspace_arguments:
                reference = addnodes.pending_xref(
                    '', refdomain='py', reftype='attr', refspecific=False,
                    reftarget=f'pyarts3.workspace.Workspace.{name}', refwarn=True,
                )
                reference += nodes.Text(f'self.{name}')
                default.children.clear()
                default += reference
            elif generic_defaults[name] is not None:
                value = generic_defaults[name].value
                default.children.clear()
                default += nodes.Text(' '.join(repr(value).splitlines()))
            else:
                # std::nullopt means required, despite the binding's None default.
                equals = next(i for i, child in enumerate(parameter.children)
                              if isinstance(child, addnodes.desc_sig_operator)
                              and child.astext() == '=')
                del parameter.children[equals:]
                if isinstance(parameter.children[-1], addnodes.desc_sig_space):
                    parameter.pop()

            # None is a binding sentinel, not an effective argument type.
            for annotation in parameter.findall(addnodes.desc_sig_name):
                children = annotation.children
                if (len(children) >= 5
                        and isinstance(children[-1], addnodes.pending_xref)
                        and children[-1].get('reftarget') == 'None'
                        and children[-3].astext() == '|'):
                    del children[-4:]


def setup(app):
    app.connect('doctree-read', format_workspace_signatures)
    return {'version': '2', 'parallel_read_safe': True, 'parallel_write_safe': True}
