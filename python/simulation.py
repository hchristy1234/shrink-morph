import shrink_morph_py
import numpy as np
import polyscope as ps
import polyscope.imgui as gui
import tkinter
from tkinter import filedialog
import triangle as tr
from svgpathtools import parse_path
import xml.etree.ElementTree as ET

root = tkinter.Tk()
root.withdraw()

lambda1 = 1
lambda2 = 0.94
thickness = 0.4

# Initialize polyscope
ps.init()
ps.set_give_focus_on_show(True)
ps.set_ground_plane_mode("shadow_only")


def callback():
  global lambda1, lambda2, P, F, theta

  if gui.Button("Select File"):
      filename = filedialog.askopenfilename(defaultextension=".obj", filetypes=[("Wavefront OBJ", "*.obj")])
      # filename = "data/circle_60.svg"

      # Parse the SVG file
      tree = ET.parse(filename)
      root = tree.getroot()

      # Handle namespaces if present
      ns = {'svg': 'http://www.w3.org/2000/svg'}

      # Extract the 'd' attribute of the first <path> element
      path_elem = root.find('.//svg:path', ns)
      if path_elem is not None:
          svg_path_data = path_elem.attrib['d']
          print(svg_path_data)
      else:
          print("No <path> element found.")
      
      # Parse the path
      path = parse_path(svg_path_data)

      # Extract vertices
      points = []
      for segment in path:
          points.append((segment.start.real, segment.start.imag))
      # Ensure it's closed
      if points[0] != points[-1]:
          points.append(points[0])

      # Convert to numpy array
      vertices = np.array(points)

      # Triangulate using the Triangle library
      A = {'vertices': vertices[:-1], 'segments': [[i, (i + 1) % (len(vertices) - 1)] for i in range(len(vertices) - 1)]}
      B = tr.triangulate(A, 'pqa0.1')

      zeros = np.zeros((len(B['vertices']), 1))
      P = np.hstack((np.array(B['vertices']), zeros))
      F = np.array(B['triangles'])
      theta = np.zeros(P.shape[0])

      # enable auto centering and scaling
      ps.set_autocenter_structures(True)
      ps.set_autoscale_structures(True)

      ps.register_surface_mesh("Mesh", P, F, edge_width=1, color=(42/255, 53/255, 213/255))
      ps.get_surface_mesh("Mesh").add_scalar_quantity("stretch orientation", theta, defined_on='vertices', enabled=True, vminmax=(-np.pi/2, np.pi/2))
        


  gui.PushItemWidth(70)
  _, lambda1 = gui.InputDouble("Lambda1", lambda1, format="%.2f")
  _, lambda2 = gui.InputDouble("Lambda2", lambda2, format="%.2f")

  if gui.Button("Simulate"):
    V = P.copy()
    shrink_morph_py.simulation(V, P[:,:2], F, theta, 10, lambda1, lambda2, thickness, 1000, 1e-6)

    ps.get_surface_mesh("Parameterization").set_enabled(False)
    ps.register_surface_mesh("Simulation", V, F, edge_width=1, color=(42/255, 53/255, 213/255))



ps.set_user_callback(callback)


ps.show()