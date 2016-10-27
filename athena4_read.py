"""
Read Athena4.2 output data files.
"""

# Python modules
import numpy as np


#=======================================================================================

def vtk(filename):
  """Read .vtk files and return dict of arrays of data."""

  # Python module
  import struct

  # Read raw data
  with open(filename, 'r') as data_file:
    raw_data = data_file.read()

  # Skip header
  current_index = 0
  current_char = raw_data[current_index]

  # Skip the first line
  # #vtk DataFile Version 3.0
  while current_char != '\n':
    current_index += 1
    current_char = raw_data[current_index]

  # Extract time info from the second line after time=...
  # CONSERVED vars at time= 1.539383e+03, level= 0, domain= 0
  while current_char != '=':
    current_index += 1
    current_char = raw_data[current_index]
  stime = ""
  while current_char != ',':
    current_index += 1
    current_char = raw_data[current_index]
    stime += current_char

  current_index += 1
  current_char = raw_data[current_index]
  time = float(stime[:-1])
  print 'time = ',time

  while current_char != '\n':
    current_index += 1
    current_char = raw_data[current_index]

  current_index += 1
  # Function for skipping though the file
  def skip_string(expected_string):
    expected_string_len = len(expected_string)
    if raw_data[current_index:current_index+expected_string_len] != expected_string:
      raise AthenaError('File not formatted as expected')
    return current_index+expected_string_len

  # Read metadata
  #BINARY
  #DATASET STRUCTURED_POINTS
  #DIMENSIONS 129 513 513
  current_index = skip_string('BINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS ')
  end_of_line_index = current_index + 1
  while raw_data[end_of_line_index] != '\n':
    end_of_line_index += 1
  face_dimensions = map(int, raw_data[current_index:end_of_line_index].split(' '))
  print 'face_dimensions = ',face_dimensions
  #cell_dimensions = [val-1 if val >1 else val for val in face_dimensions]
  #print 'cell_dimensions = ',cell_dimensions

  current_index = end_of_line_index + 1
  # Read interface locations
  #ORIGIN -5.000000e-01 -2.000000e+00 -2.000000e+00
  current_index = skip_string('ORIGIN')+1
  end_of_line_index = current_index + 1
  while raw_data[end_of_line_index] != '\n':
    end_of_line_index += 1
  print raw_data[current_index:end_of_line_index]
  box_origin = map(float, raw_data[current_index:end_of_line_index].split(' '))
  print 'box_origin = ', box_origin

  current_index = end_of_line_index + 1
  #SPACING 7.812500e-03 7.812500e-03 7.812500e-03
  current_index = skip_string('SPACING')+1
  end_of_line_index = current_index + 1
  while raw_data[end_of_line_index] != '\n':
    end_of_line_index += 1
  print raw_data[current_index:end_of_line_index]
  grid_spacing = map(float, raw_data[current_index:end_of_line_index].split(' '))
  print 'grid_spacing = ', grid_spacing


  current_index = end_of_line_index + 1
  # Prepare to read quantities defined on grid
  cell_dimensions = np.array([max(dim-1,1)
      for dim in face_dimensions])
  #CELL_DATA 33554432
  num_cells = cell_dimensions.prod()
  current_index = skip_string('CELL_DATA {0}\n'.format(num_cells))
  end_of_line_index = current_index + 1
  print 'cell_dimensions = ',cell_dimensions[::-1],' num_cells = ',num_cells


  # construct the cell centered grid based on origin, spacing and cell_dimensions
  x = box_origin[0] + np.arange(cell_dimensions[0])*grid_spacing[0] + 0.5*grid_spacing[0]
  y = box_origin[1] + np.arange(cell_dimensions[1])*grid_spacing[1] + 0.5*grid_spacing[1]
  z = box_origin[2] + np.arange(cell_dimensions[2])*grid_spacing[2] + 0.5*grid_spacing[2]
  #print "x-grid: "
  #print x
  #print "y-grid: "
  #print y
  #print "z-grid: "
  #print z


  #if raw_data[current_index:current_index+1] == '\n':
  #  current_index = skip_string('\n')  # extra newline inserted by join script

  data = {}

  # Function for reading scalar data
  def read_cell_scalars():
    begin_index = skip_string('SCALARS ')
    end_of_word_index = begin_index + 1
    while raw_data[end_of_word_index] != ' ':
      end_of_word_index += 1
    array_name = raw_data[begin_index:end_of_word_index]
    print 'loading array_name = ',array_name
    string_to_skip = 'SCALARS {0} float\nLOOKUP_TABLE default\n'.format(array_name)
    begin_index = skip_string(string_to_skip)
    format_string = '>' + 'f'*num_cells
    end_index = begin_index + struct.calcsize('f')*num_cells
    print 'is,ie = ',begin_index, end_index
    data[array_name] = struct.unpack(format_string, raw_data[begin_index:end_index])
    dimensions = tuple(cell_dimensions[::-1])
    data[array_name] = np.array(data[array_name]).reshape(dimensions)
    #dimensions = tuple(cell_dimensions)
    #data[array_name] = np.array(data[array_name]).reshape(dimensions)
    print 'data['+array_name+'].shape= ',data[array_name].shape
    return end_index #+1

  # Function for reading vector data
  def read_cell_vectors():
    begin_index = skip_string('VECTORS ')
    end_of_word_index = begin_index + 1
    while raw_data[end_of_word_index] != ' ':
      end_of_word_index += 1
    array_name = raw_data[begin_index:end_of_word_index]
    print 'loading array_name = ',array_name
    string_to_skip = 'VECTORS {0} float\n'.format(array_name)
    ###array_name = array_name[:-6]  # remove ' float'
    begin_index = skip_string(string_to_skip)
    format_string = '>' + 'f'*num_cells*3
    end_index = begin_index + struct.calcsize('f')*num_cells*3
    data[array_name] = struct.unpack(format_string, raw_data[begin_index:end_index])
    dimensions = tuple(np.append(cell_dimensions[::-1],3))
    data[array_name] = np.array(data[array_name]).reshape(dimensions)
    #dimensions = tuple(np.append(3,cell_dimensions))
    data[array_name] = np.array(data[array_name]).reshape(dimensions)
    return end_index #+1

  # Read quantities defined on grid
  #SCALARS density float
  #LOOKUP_TABLE default
  #VECTORS momentum float
  #VECTORS cell_centered_B float
  while current_index < len(raw_data):
    expected_string = 'SCALARS'
    expected_string_len = len(expected_string)
    print raw_data[current_index:current_index+expected_string_len]
    if raw_data[current_index:current_index+expected_string_len] == expected_string:
      print "start loading scalars !! "
      current_index = read_cell_scalars()
      continue

    expected_string = 'VECTORS'
    expected_string_len = len(expected_string)
    print raw_data[current_index:current_index+expected_string_len]
    if raw_data[current_index:current_index+expected_string_len] == expected_string:
      print "start loading vectors !! "
      current_index = read_cell_vectors()
      continue
#  raise Athena4Error('File not formatted as expected')



  return time,x,y,z,data

#=======================================================================================

class Athena4Error(RuntimeError):
  """General exception class for Athena4 read functions."""
  pass
