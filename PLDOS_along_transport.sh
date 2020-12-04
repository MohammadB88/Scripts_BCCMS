#!/bin/usr/

################################################################################################
# DOS projected on atoms from leyers of the left electrode to the right electrode, TO plot PDOS along the transport direction:

########################################################
# read PLDOS of layers before V_Mo on the left electrode
for (( l=4; l<9; ++l))

do 
        echo $l
        #echo $((24*l))
        beg_atom=$(($((24*$((l-1))))+1))
        end_atom=$((24*l))
        echo ${beg_atom} ${end_atom}
        sdata MoS2_dev.TBT.nc --atom ${beg_atom}-${end_atom} --dos --ados Left --out PDOS_along_Z/dos_${l}L_${beg_atom}-${end_atom}.dat

done

###############################################################
# read PLDOS of the layer containing V_Mo on the left electrode
l_d_left=9
echo 'This layer has vacancy at the LEFT electrode!'
echo $l_d_left
beg_atom=$(($((24*$((l_d_left -1))))+1))
end_atom=$(($((24*l_d_left))-1)) # -1 weil one V_Mo
echo ${beg_atom} ${end_atom}
sdata MoS2_dev.TBT.nc --atom ${beg_atom}-${end_atom} --dos --ados Left --out PDOS_along_Z/dos_${l_d_left}L_${beg_atom}-${end_atom}.dat

#######################################################
# read PLDOS of layers after V_Mo on the left electrode
# and layers before V_Mo on the right electrode
for (( l=l_d_left+1; l<29; ++l))

do 
        echo $l
        #echo $((24*l))
        beg_atom=$(($(($((24*$((l-1))))+1))-1)) # -1 weil one V_Mo
        end_atom=$(($((24*l))-1)) # -1 weil one V_Mo
        echo ${beg_atom} ${end_atom}
        sdata MoS2_dev.TBT.nc --atom ${beg_atom}-${end_atom} --dos --ados Left --out PDOS_along_Z/dos_${l}L_${beg_atom}-${end_atom}.dat

done

################################################################
# read PLDOS of the layer containing V_Mo on the right electrode
l_d_right=29
echo 'This layer has vacancy at the RIGHT electrode!'
echo $l_d_right
beg_atom=$(($(($((24*$((l_d_right -1))))+1))-1)) # -1 weil one V_Mo
end_atom=$(($((24*l_d_right))-2)) # -2 weil one V_Mo + one V_Mo
echo ${beg_atom} ${end_atom}
sdata MoS2_dev.TBT.nc --atom ${beg_atom}-${end_atom} --dos --ados Left --out PDOS_along_Z/dos_${l_d_right}L_${beg_atom}-${end_atom}.dat

############################################################
# read PLDOS of the layers after V_Mo on the right electrode
for (( l=l_d_right+1; l<35; ++l))

do 
        echo $l
        #echo $((24*l))
        beg_atom=$(($(($((24*$((l-1))))+1))-2)) # -2 weil one V_Mo + one V_Mo
        end_atom=$(($((24*l))-2)) # -2 weil one V_Mo + one + V_Mo
        echo ${beg_atom} ${end_atom}
        sdata MoS2_dev.TBT.nc --atom ${beg_atom}-${end_atom} --dos --ados Left --out PDOS_along_Z/dos_${l}L_${beg_atom}-${end_atom}.dat

done
