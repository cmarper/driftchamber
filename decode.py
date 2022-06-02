#!/usr/bin/env python

"""
Script to convert binary format to numpy arrays for data acquired with DRS4 evaluation boards.
http://www.psi.ch/drs/evaluation-board
http://www.psi.ch/drs/DocumentationEN/manual_rev40.pdf
By: Jonas Rembser, Gregor Kasieczka, Dmitry Hits, Cristina Martin Perez
"""

from sys import argv, exit
from numpy import array, uint32, cumsum, roll, zeros, float32, arange
from struct import unpack
import matplotlib.pyplot as plt

# To run: python decode.py filename.dat

# Inputs
if not len(argv) == 2:
    print("Wrong number of arguments!")
    print("Usage: python decode.py filename.dat")
    print("Exiting...")
    exit()

input_filename = argv[1]
f = open(input_filename, "rb")

# Board and event serial number
board_serials = []
event_serial = array([0], dtype=uint32)

# Total number of channels and boards
n_ch = 0
n_boards = 0

# Time and voltage numpy arrays (one per channel)
channels_t = []
channels_v = []

# List of numpy arrays to store the events and time bin information
events = []
timebins = []

# Extract time information for each DRS4 cell
while True:
    header = f.read(4)
    # Skip the initial time header
    if header == b"TIME":
        continue
    elif header.startswith(b"C"):
        n_ch = n_ch + 1
        # Create empty arrays for voltage and time
        channels_t.append(zeros(1024, dtype=float32))
        channels_v.append(zeros(1024, dtype=float32))
        # Write timebins to numpy array
        timebins.append(array(unpack('f'*1024, f.read(4*1024))))

    # Increment the number of boards when seeing a new serial number
    # and store the serial numbers in the board serial numbers vector
    elif header.startswith(b"B#"):
        board_serial = unpack(b'H', header[2:])[0]
        board_serials.append(board_serial)
        n_boards = n_boards + 1

    # End the loop if header is not CXX or a serial number
    elif header == b"EHDR":
        break

# Main loop
current_board = 0
tcell = 0 # current trigger cell
t_00 = 0 # time in first cell in first channel (for alignment)
is_new_event = True

print('\n')
info_string = "Reading events measurend with {0} channels on {1} board(s)..."
print(info_string.format(n_ch, n_boards))

while True:
    # Start of Event
    if is_new_event:
        event_serial[0] = unpack("I", f.read(4))[0]
        if event_serial[0]%10 == 0:
            print("Running on event : ", event_serial[0])
        is_new_event = False
        # Set the timestamp, where the milliseconds need to be converted to
        # nanoseconds to fit the function arguments
        dt_list = unpack("H"*8, f.read(16))

        # Fluff the serial number and read in trigger cell
        fluff = f.read(4)
        tcell = unpack('H', f.read(4)[2:])[0]
        # Reset current board number
        current_board = 0
        continue

    # Read the header, this is either
    #  EHDR -> finish event
    #  C00x -> read the data
    #  ""   -> end of file
    header = f.read(4)

    # Handle next board
    if header.startswith(b"B#"):
        current_board = current_board + 1
        tcell = unpack(b'H', f.read(4)[2:])[0]
        continue

    # End of Event
    elif header == b"EHDR":
        events.append(array(list(zip(channels_t[1],channels_v[1]))))
        is_new_event = True

    # Read and store data from a channel read from the header
    elif header.startswith(b'C'):
        chn_i = int(header.decode('ascii')[-1]) + current_board * 4 # channel number - 1
        scaler = unpack('I', f.read(4))
        voltage_ints = unpack(b'H'*1024, f.read(2*1024)) # the voltage info is 1024 floats with 2-byte precision

        # Calculate precise timing using the time bins and trigger cell
        # Sum up the times of all cells starting from the trigger cell to the i_th cell 
        # and select only even members, as the amplitude of the adjacent cells are averaged.
        timebins_full = list(roll(timebins[chn_i-1], -tcell))+list(roll(timebins[chn_i-1], -tcell))
        t = cumsum(timebins_full)[::2]
        # Find the time of the first cell for alignment correction
        t_0 = t[(1024-tcell)%1024]
        if chn_i % 4 == 1:
            t_00 = t_0
        # Align all channels with the first channel
        t = t - (t_0 - t_00) # correction

        # Output list of numpy arrays (to be plotted or used in the further analysis).
        # Time vector per event: channels_t[channelnumber][cellnumber]
        # Voltage vector per event: channels_v[channelnumber][cellnumber]
        # NOTE: channels 0-3 are for the first board,  channels 4-7 are for the second board
        for i, x in enumerate(voltage_ints): # i is the cell number, x is the voltage value of the sample in ADC channels
            channels_v[chn_i-1][i] = ((x / 65535.) - 0.5)
            channels_t[chn_i-1][i] = t[i]

        # Print times and voltages per event per channel
        #print('  Channel: ',chn_i-1)
        #print('    Times: ', channels_t[chn_i-1])
        #print('    Voltages: ', channels_v[chn_i-1])

    # End of File
    elif header == b"":
        break

f.close()
print('\n')
