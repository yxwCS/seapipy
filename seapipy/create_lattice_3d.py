import numpy as np
from scipy.spatial import Voronoi, ConvexHull
from scipy.spatial import distance


class create_lattice_3d:
    def __init__(self, number_cells_x, number_cells_y, number_cells_z, standard_deviation=0.15, spatial_step=20):
        self.number_cells_x = number_cells_x
        self.number_cells_y = number_cells_y
        self.number_cells_z = number_cells_z
        self.standard_deviation = standard_deviation
        self.spatial_step = spatial_step

    def generate_cube_seeds(self):
        """
        Generate seeds for tessellation in a cubic grid.
        """
        grid_values = [
            [
                (i + np.random.normal(0, self.standard_deviation)) * self.spatial_step,
                (j + np.random.normal(0, self.standard_deviation)) * self.spatial_step,
                (k + np.random.normal(0, self.standard_deviation)) * self.spatial_step
            ]
            for k in range(self.number_cells_z) for j in range(self.number_cells_y) for i in range(self.number_cells_x)
        ]
        return grid_values

    

    def generate_edges_and_faces_from_vertices_v2(self, faces_dict):
        edges = {}
        edge_map = {}
        edge_num = 1
        edge_faces = {}

        for index, face in faces_dict.items():
            edge_face = []
            for ii in range(len(face)):
                v0 = face[ii]
                v1 = face[(ii + 1) % len(face)]
                edge = (v0, v1) if v0 < v1 else (v1, v0)
                if edge not in edge_map:
                    edges[edge_num] = edge
                    edge_map[edge] = edge_num
                    edge_map[(edge[1], edge[0])] = -edge_num
                    edge_num += 1
                edge_face.append(edge_map[(v0, v1)])
            edge_faces[index] = edge_face
        return edges, edge_faces

    def lloyd_relax(self, points, steps=1):
        for _ in range(steps):
            vor = Voronoi(points)
            new_points = []
            for region in vor.regions:
                if not region:
                    continue
                if -1 in region:
                    continue
                vertices = vor.vertices[region]
                centroid = np.mean(vertices, axis=0)  # calculate the centroid of each polygon
                new_points.append(centroid)
            if len(new_points) >= 5:  # Check if enough points exist
                points = np.array(new_points)  # Only update if there are enough points
            else:
                break  # Exit the loop if not enough points
        return points

    def find_distant_faces(self, faces, vertices, max_distance=100) -> dict:
        distant_faces_index = []

        for index, face in faces.items():
            if -1 in face:
                continue

            is_distant = False

            for i in range(len(face)):
                vertex_index1 = face[i]
                vertex_index2 = face[(i + 1) % len(face)]
                point1 = vertices[vertex_index1]
                point2 = vertices[vertex_index2]
                dis = distance.euclidean(point1, point2)

                if dis > max_distance:
                    is_distant = True
            
                    break

            if is_distant:
                distant_faces_index.append(index)
        return distant_faces_index

    def bodies_generator(self, faces_index_map, bodies, invalid_faces_index=[]):
        bodies_by_faces = []
        body_vertices_set = set()
        for body in bodies:
            if body:
                body_vertices_set = set(body)
                body_temp = []
                for index, face in faces_index_map.items():
                    face_set = set(face)
                    if face_set.issubset(body_vertices_set):
                        body_temp.append(index)
                bodies_by_faces.append(body_temp)

        final_bodies = []
        for body in bodies_by_faces:
            if body:
                if not any(face in invalid_faces_index for face in body):
                    final_bodies.append(body)
        return final_bodies

    def filter_faces_by_bodies(self, filtered_faces, bodies_by_faces):
        final_faces = {}
        for face, i in filtered_faces.items():
            
            for body in bodies_by_faces:
                if i in body:
                    
                    final_faces[i] = face
        return final_faces


    def remove_invalid_faces(self, faces_index, invalid_faces):
        valid_faces = {vertices: number for number, vertices in faces_index.items() if number not in invalid_faces}
        return valid_faces


    def filter_vertices(self, all_vertices_dict, edges):
        ref_vertices = set()
        for v in edges.values():
            ref_vertices.update(v)

        
        valid_vertices_dict = {vertex: data for vertex, data in all_vertices_dict.items() if vertex in ref_vertices}

        return valid_vertices_dict


    def rearrange_faces_bodies(self, vertices, faces, bodies):
        vertex_dict = {i: vertices[i] for i in range(len(vertices))}

        used_vertices = set()
        for face in faces:
            used_vertices.update(face)

        for body in bodies:
            used_vertices.update(body)

        new_index = {}
        new_vertices = []
        for old_idx in used_vertices:
            new_index[old_idx] = len(new_vertices) + 1
            new_vertices.append(vertex_dict[old_idx - 1])

        new_faces = [[new_index[idx] for idx in face] for face in faces]
        new_bodies = [[new_index[idx] for idx in body] for body in bodies]

        return new_vertices, new_faces, new_bodies

    def volume_cal(self, vertices):
        vertices = np.array(vertices)
        hull = ConvexHull(vertices)
        volume = hull.volume
        return volume

    def bodies_volume_assign(self, bodies_by_faces, vertices_index_map, faces_by_vertices):
        volume_list = []

        for body in bodies_by_faces:
            vertices = set()
            for face in body:
                vertices_of_face = faces_by_vertices[face]
                vertices.update(vertices_of_face)

            vertices_cor = [vertices_index_map[vertex] for vertex in vertices]
            print("Vertices cor", vertices_cor)
            volume = self.volume_cal(vertices_cor)
            volume_list.append(volume)

        return volume_list

    def rearrange_vertex_indices(self, valid_vertices, edges, faces, bodies):
        sorted_valid_vertices = sorted(valid_vertices.keys())
        old_to_new_v = {old_index: new_index for new_index, old_index in enumerate(sorted_valid_vertices, start=1)}

        new_vertices = {old_to_new_v[old]: valid_vertices[old] for old in sorted_valid_vertices}

        new_edges = {edge_id: (old_to_new_v[edge[0]], old_to_new_v[edge[1]]) for edge_id, edge in edges.items()}

        old_to_new_f = {face_id: new_face_id for new_face_id, face_id in enumerate(faces.keys(), start=1)}
        new_faces = {old_to_new_f[face_id]: face for face_id, face in faces.items()}

        new_bodies = [[old_to_new_f[face_id] for face_id in body if face_id in old_to_new_f] for body in bodies]

        return new_vertices, new_edges, new_faces, new_bodies


    def calculate_normal(self, v1, v2, v3):
        vector1 = np.array(v2) - np.array(v1)
        vector2 = np.array(v3) - np.array(v1)
        normal = np.cross(vector1, vector2)
        return normal

    def calculate_centroid(self, vertices):
        centroid = np.mean(vertices, axis=0)
        return centroid

    def orient_body_faces(self, body, vertices_index_map, edges_by_vertices, faces_by_edges):
        body_vertices = set()
        for face in body:
            for edge in faces_by_edges[abs(face)]:
                body_vertices.update(edges_by_vertices[abs(edge)])
        body_vertices = list(body_vertices)
        body_centroid = self.calculate_centroid([vertices_index_map[vertex] for vertex in body_vertices])

        oriented_body = []
        for face in body:
            face_edges = faces_by_edges[abs(face)]
            face_vertices = set()
            for edge in face_edges:
                face_vertices.update(edges_by_vertices[abs(edge)])
            face_vertices = list(face_vertices)
            v1, v2, v3 = (vertices_index_map[face_vertices[i]] for i in range(3))

            normal = self.calculate_normal(v1, v2, v3)
            centroid_vector = np.array(body_centroid) - np.array(v1)
            dot_product = np.dot(normal, centroid_vector)

            if dot_product > 0:  # inward normal
                oriented_body.append(-face)
            else:
                oriented_body.append(face)
        return oriented_body

    def orient_bodies(self, bodies_by_faces, vertices_index_map, edges_by_vertices, faces_by_edges):
        oriented_bodies = []
        for body in bodies_by_faces:
            oriented_body = self.orient_body_faces(body, vertices_index_map, edges_by_vertices, faces_by_edges)
            oriented_bodies.append(oriented_body)
        return oriented_bodies

    def create_lattice(self):
        seeds = self.generate_cube_seeds()
        vor = Voronoi(seeds)
        vertices = vor.vertices
        vertices_index_map = {i+1: tuple(vertex) for i, vertex in enumerate(vor.vertices)}
        pos_faces = [[vertex + 1 for vertex in face] for face in vor.ridge_vertices if -1 not in face]
        all_faces_index = {i + 1: tuple(face) for i, face in enumerate(pos_faces)}
        bodies = [[vertex + 1 for vertex in region] for region in vor.regions if -1 not in region]
        invalid_faces = self.find_distant_faces(all_faces_index, vertices_index_map, max_distance=40)
        bodies_by_faces = self.bodies_generator(all_faces_index, bodies, invalid_faces)

        filtered_faces = self.remove_invalid_faces(all_faces_index, invalid_faces)

        final_faces = self.filter_faces_by_bodies(filtered_faces, bodies_by_faces)

        edges, edge_faces = self.generate_edges_and_faces_from_vertices_v2(final_faces)

        ref_vertices = self.filter_vertices(vertices_index_map, edges)

        volume_list = self.bodies_volume_assign(bodies_by_faces, vertices_index_map, all_faces_index)

        new_vertices, new_edges, new_faces, new_bodies = self.rearrange_vertex_indices(ref_vertices, edges, edge_faces, bodies_by_faces)
        oriented_bodies = self.orient_bodies(new_bodies, new_vertices, new_edges, new_faces)
        return new_vertices, new_faces, oriented_bodies, volume_list

