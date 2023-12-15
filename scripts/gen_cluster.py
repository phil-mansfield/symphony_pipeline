names = ["Halo156", "Halo175", "Halo200", "Halo211", "Halo213", "Halo222", "Halo225", "Halo266", "Halo274", "Halo277", "Halo282", "Halo293", "Halo304", "Halo305", "Halo306", "Halo308", "Halo317", "Halo321", "Halo324", "Halo326", "Halo335", "Halo337", "Halo339", "Halo345", "Halo346", "Halo348", "Halo349", "Halo352", "Halo354", "Halo358", "Halo360", "Halo361", "Halo366", "Halo367", "Halo377", "Halo378", "Halo385", "Halo386", "Halo387", "Halo390", "Halo391", "Halo394", "Halo400", "Halo407", "Halo409", "Halo415", "Halo416", "Halo419", "Halo428", "Halo429", "Halo436", "Halo437", "Halo441", "Halo445", "Halo447", "Halo448", "Halo452", "Halo454", "Halo455", "Halo456", "Halo461", "Halo462", "Halo465", "Halo471", "Halo472", "Halo474", "Halo475", "Halo476", "Halo478", "Halo479", "Halo480", "Halo483", "Halo489", "Halo494", "Halo502", "Halo517", "Halo518", "Halo522", "Halo529", "Halo544", "Halo545", "Halo546", "Halo561", "Halo572", "Halo574", "Halo595", "Halo600", "Halo604", "Halo629", "Halo631", "Halo639", "Halo645", "Halo653", "Halo734", "Halo372", "Halo425"]

def r(l, h): return list(range(l, h+1))

idxs = r(0, 12) + r(14, 45) + r(48, 48) + r(51, 55) + r(57, 60) + r(64, 75) + r(77, 77) + r(79, 84) + r(86, 86) + r(88, 88) + r(90, 91) + r(93, 95)
print(idxs)
print("[", end="")
for i in idxs:
    print("'%s'" % names[i], end=", ")
print("]")
