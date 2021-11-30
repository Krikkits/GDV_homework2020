#----------------------------------------------------------
# File export_ra2.py
# Simple ra2 exporter which dumps 3 float values for each of the 3 vertices of each triangle in a row.
#----------------------------------------------------------
import bpy, os
import struct
 
def export_ra2(filepath, ob, rot90, scale):
    name = os.path.basename(filepath)
    realpath = os.path.realpath(os.path.expanduser(filepath))
    meshfile = open(realpath, 'wb')    
    print('Exporting %s' % realpath)
 
    if not ob or ob.type != 'MESH':
        raise NameError('Cannot export: active object %s is not a mesh.' % ob)
    me = ob.data
    
    for f in me.polygons:
        vs = f.vertices
        vertex1 = scale * me.vertices[f.vertices[0]].co
        vertex2 = scale * me.vertices[f.vertices[1]].co
        vertex3 = scale * me.vertices[f.vertices[2]].co
        if rot90:
            meshfile.write(struct.pack('fff', vertex1[0], vertex1[2], -vertex1[1]))
            meshfile.write(struct.pack('fff', vertex2[0], vertex2[2], -vertex2[1]))
            meshfile.write(struct.pack('fff', vertex3[0], vertex3[2], -vertex3[1]))
        else:
            meshfile.write(struct.pack('fff', vertex1[0], vertex1[1], vertex1[2]))
            meshfile.write(struct.pack('fff', vertex2[0], vertex2[1], vertex2[2]))
            meshfile.write(struct.pack('fff', vertex3[0], vertex3[1], vertex3[2]))
        if len(f.vertices) == 4:
            vertex4 = scale * me.vertices[f.vertices[3]].co
            if rot90:
                meshfile.write(struct.pack('fff', vertex1[0], vertex1[2], -vertex1[1]))
                meshfile.write(struct.pack('fff', vertex3[0], vertex3[2], -vertex3[1]))
                meshfile.write(struct.pack('fff', vertex4[0], vertex4[2], -vertex4[1]))
            else:
                meshfile.write(struct.pack('fff', vertex1[0], vertex1[1], vertex1[2]))
                meshfile.write(struct.pack('fff', vertex3[0], vertex3[1], vertex3[2]))
                meshfile.write(struct.pack('fff', vertex4[0], vertex4[1], vertex4[2]))
 
    print('%s successfully exported' % realpath)
    meshfile.close()
    return
