SRC = mfiu.sv trip_row.sv trgt_sub.sv trgs_row.sv trapezoid_top.sv
all: unit sweep
unit: unit-dense unit-trip unit-trgt unit-trgs
unit-dense:
	iverilog -g2012 -o tb_dense.vvp mfiu.sv trip_row.sv tb_dense.sv && vvp tb_dense.vvp | tail -3
unit-trip:
	iverilog -g2012 -o tb_trip.vvp mfiu.sv trip_row.sv tb_trip.sv && vvp tb_trip.vvp | tail -3
unit-trgt:
	iverilog -g2012 -o tb_trgt.vvp trgt_sub.sv tb_trgt.sv && vvp tb_trgt.vvp | tail -3
unit-trgs:
	iverilog -g2012 -o tb_trgs.vvp trgs_row.sv tb_trgs.sv && vvp tb_trgs.vvp | tail -3
sweep:
	iverilog -g2012 -o tb_top.vvp $(SRC) tb_top.sv && vvp tb_top.vvp | tail -3

