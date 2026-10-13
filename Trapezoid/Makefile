SRC = mfiu.sv trip_row.sv trgt_sub.sv trgs_row.sv trapezoid_top.sv
all: unit sweep
unit:
	iverilog -g2012 -o tb_trip.vvp mfiu.sv trip_row.sv tb_trip.sv && vvp tb_trip.vvp | tail -3
sweep:
	iverilog -g2012 -o tb_top.vvp $(SRC) tb_top.sv && vvp tb_top.vvp | tail -3
