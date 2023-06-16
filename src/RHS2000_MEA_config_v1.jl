"""
        this function takes a string of channel sequences and matches the right
        sequence of channels from a 128 channel Intan data aquistion board (model: RHS2000)
        according to the multielectrode array (MEA) configuration.

        The string (probes) contains a sequence of 8 characters, depending on which port were used
        in the experiment. Each port carries 32 channels of data:

        00 = port not connected to a MEA
        s1 = surface array (32 channels)
        p1 = penetrating electrode of 16 x 8
        p4 = 4 shank of 32 channels

        eg: "p4p4p4p4" = a 4 shank probe connected to all 4 ports
        "s100p1p1" = surface probe in port 1, nothing in port 2, and single shank 32
        channel probes in ports 3 and 4.

        possible combos so far:

        s1s1p1p1
        s1p10000
        p4p4p4p4

        note that mixing p1 with p4 is problematic and will require custom configuration
        because a single shank of p4 is not mapped across a single port.
"""
function rhsMEA_config(probes)

        chConfig=zeros(16,8)
        # port_off=zeros(16,2)

        s1base1 = # this electrode config makes no sense, so will use p1 config for s1 until resolved
        [9	24	21	25	8
        15	17	28	32	2
        10	23	16	26	7
        14	18	1	31	3
        11	22	12	27	6
        13	19	1	30	4
        1	29	5	20	1]

        p1base1 =
        [16  32
        21  26
        18  31
        13   7
        15   3
        20   4
        10  30
        14   5
        19  29
        23  28
        12   2
        17   1
        22   6
        11  27
        24   8
        9  25]

        p4_all=
        [120 63 33 106 87 32 2 73;
        118 61 35 108 85 30 4 75;
        116 59 37 110 83 28 6 77;
        114 57 39 112 81 26 8 79;
        121 50 48 103 90 17 15 72;
        123 52 46 101 92 19 13 70;
        125 54 44 99 94 21 11 68;
        127 56 42 97 96 23 9 66;
        128 55 41 98 95 24 10 65;
        126 53 43 100 93 22 12 67;
        124 51 45 102 91 20 14 69;
        122 49 47 104 89 18 16 71;
        119 64 34 105 88 31 1 74;
        117 62 36 107 86 29 3 76;
        115 60 38 109 84 27 5 78;
        113 58 40 111 82 25 7 80]

        if probes[1:2] == "p4" # port A connected
                chConfig[:, 7] = p4_all[:,7]
                chConfig[:, 6] = p4_all[:,6]
        # elseif  probes[1:2] == "p1" || probes[1:2] == "s1" # port A connected
        elseif probes[1:2] == "00" # port A connected
                chConfig[:, 1:2] = zeros(16,2)
        else        chConfig[:, 1:2] = p1base1
        end

        if probes[3:4] == "p4" # port B connected
                chConfig[:, 2] = p4_all[:,2]
                chConfig[:, 3] = p4_all[:,3]
                # elseif  probes[3:4] == "p1" || probes[3:4] == "s1" # port B connected

        elseif probes[3:4] == "00" # port A connected
                chConfig[:, 3:4] = zeros(16,2)        
        else    chConfig[:, 3:4] = p1base1 .+ 32
        end

        if probes[5:6] == "p4" # port C connected
                chConfig[:, 5] = p4_all[:,5]
                chConfig[:, 8] = p4_all[:,8]
                # elseif  probes[5:6] == "p1" || probes[5:6] == "s1" # port C connected
                elseif probes[5:6] == "00" # port A connected
                chConfig[:, 5:6] = zeros(16,2)          
        else        chConfig[:, 5:6] = p1base1 .+ 64
        end

        if probes[7:8] == "p4" # port D connected
                chConfig[:, 1] = p4_all[:,1]
                chConfig[:, 4] = p4_all[:,4]
                # elseif  probes[7:8] == "p1" || probes[7:8] == "s1" # port D connected
        elseif probes[7:8] == "00" # port A connected
                chConfig[:, 7:8] = zeros(16,2) 
        else        chConfig[:, 7:8] = p1base1 .+ 96
        end

        chConfig = round.(Int, chConfig)
        return chConfig
end
