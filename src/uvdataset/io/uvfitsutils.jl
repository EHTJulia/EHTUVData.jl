const uvfits_polid2name = Dict(
    "+1" => "I",
    "+2" => "Q",
    "+3" => "U",
    "+4" => "V",
    "-1" => "RR",
    "-2" => "LL",
    "-3" => "RL",
    "-4" => "LR",
    "-5" => "XX",
    "-6" => "YY",
    "-7" => "XY",
    "-8" => "YX",
)

const uvfits_polname2id = Dict(
    "I" => +1,
    "Q" => +2,
    "U" => +3,
    "V" => +4,
    "RR" => -1,
    "LL" => -2,
    "RL" => -3,
    "LR" => -4,
    "XX" => -5,
    "YY" => -6,
    "XY" => -7,
    "YX" => -8,
)