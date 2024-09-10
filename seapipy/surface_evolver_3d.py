
class SurfaceEvolver:
    def __init__(self, vertices, edges, faces, cells, density_values, volume_values, polygonal=True):
        self.vertices = vertices
        self.edges = edges
        self.faces = faces
        self.cells = cells

        self.density_values = {tuple(map(tuple, np.round(np.array(key), 3))): round(value, 3) for key, value in density_values.items()}
        self.volume_values = volume_values
        self.polygonal = polygonal
        self.fe_file = io.StringIO()

    def generate_fe_file(self):
        self.fe_file.write("SPACE_DIMENSION 3 \n")
        self.fe_file.write("SCALE 0.005 FIXED\n")
        # self.fe_file.write("STRING \n")
        self.fe_file.write("\n")

        # Write vertices
        self.fe_file.write("vertices \n")
        for k, v in self.vertices.items():
            self.fe_file.write(f"{k} {v[0]} {v[1]} {v[2]} \n")
        self.fe_file.write("\n")

        # Write edges
        self.fe_file.write("edges \n")
        for k, v in self.edges.items():
            self.fe_file.write(f"{k} {v[0]} {v[1]} \n")
        self.fe_file.write("\n")

        # Write faces
        self.fe_file.write("faces \n")
        for i, face in self.faces.items():  # faces variable is a dictionary
            self.fe_file.write(f"{i } {' '.join(map(str, face))} \n")
        self.fe_file.write("\n")

        # Write bodies
        self.fe_file.write("bodies \n")
        i = 1
        for k, v in enumerate(self.cells):
            if v:
                volume = self.volume_values[k]
                str_value = " ".join(str(vv) for vv in v)
                self.fe_file.write(f"{abs(i)} {str_value} volume {volume}\n")
                i += 1
        self.fe_file.write("\n")

        # Additional commands
        self.fe_file.write("read \n \n")
        self.fe_file.write("show_all_edges off \n")
        self.fe_file.write("metric_conversion off \n")
        self.fe_file.write("autorecalc on \n")
        self.fe_file.write("gv_binary off \n")
        self.fe_file.write("gravity off \n")
        self.fe_file.write("ii := 0; \n")
        self.fe_file.write("gogo := { g; r; g 1000; } \n")
        self.fe_file.write("s \n")

        return self.fe_file

    def save_fe_file(self, file_name: str):
        with open(f'{file_name}', mode='w') as f:
            print(self.fe_file.getvalue(), file=f)
        return True
