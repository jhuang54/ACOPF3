Maui2022dm_v4_v33.raw:
PSSE raw file given by Xinyang

Maui2022dm_rd_v33.raw:
Delete isolated buses in Maui2022dm_v4_v33.raw

The following three files revise Maui2022dm_rd_v33.raw, including 
(1) Maui2022dm_rd_v33_shunt.raw
(2) Maui2022dm_rd_v33_shunt_NoPhaseShift.raw
(3) Maui2022dm_rd_v33_shunt_NoPhaseShiftNoTapRatio.raw
The modified details include:
	Maui2022dm_rd_v33_shunt.raw
	convert switched shunts to generator only generating reactive power in Maui2022dm_rd_v33.raw

	Maui2022dm_rd_v33_shunt_NoPhaseShift.raw:
	(1) convert switched shunts to generator only generating reactive power
	(2) transformer angle shift:0
	in Maui2022dm_rd_v33.raw

	Maui2022dm_rd_v33_shunt_NoPhaseShiftNoTapRatio.raw:
	(1) convert switched shunts to generator only generating reactive power
	(2) transformer angle shift:0
	(3) transformer tap ratio:1


