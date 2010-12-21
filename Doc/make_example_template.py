import sys

template = """
##############################################

.. literalinclude:: ../../../%s
   :language: python
   :linenos:

"""

script_name = sys.argv[1]
file(script_name+'.rst', 'w').write(template % script_name)
