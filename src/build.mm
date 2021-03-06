pckg := "MapleCodeGenerator";
pckgname := "Maple Code Generator";

NMPCT := cat(pckg,".maple");
PackageTools:-Create(NMPCT);
PackageTools:-SetProperty(NMPCT, "Author", "Behzad Samadi");
PackageTools:-SetProperty(NMPCT, "Item List", "true");
PackageTools:-SetProperty(NMPCT, "X-CloudGroup", "packages"); 
PackageTools:-SetProperty(NMPCT, "X-CloudId", "6323129933627392"); 
PackageTools:-SetProperty(NMPCT, "X-CloudURL", "https://maple.cloud");
PackageTools:-SetProperty(NMPCT, "X-CloudXId", "behzad.samadi@gmail.com");

# Overview
PackageTools:-AddAttachment(NMPCT, "/Overview.mw" = "../doc/overview.mw");
PackageTools:-AddAttachment(NMPCT, "/Examples/NMPC.mw" = "../doc/examples-nmpc.mw");

read cat(pckg, ".mpl");
savelib(convert(pckg,name),NMPCT);