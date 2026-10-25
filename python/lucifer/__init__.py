"""The Python side of Lucifer.

deck    writes the TOML deck the program reads, from flat key = value text or a mapping.
schema  loads the JSON Schema the program writes with lucifer --schema.
io      reads the stats file and the openPMD beam and field dumps.

Nothing here holds a default the program does not, generates anything but the TOML the
program reads, or checks what the program refuses itself. The program is the source of
truth, and this package makes decks and reads results.
"""
