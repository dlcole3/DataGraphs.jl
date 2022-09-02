using Colors, TestImages, Images


img = Images.load("./examples/Bucky_Badger.jpg")
#picture from https://www.wikiwand.com/en/Bucky_Badger

imgg = Gray.(img)

mat = convert(Array{Float64}, imgg)
