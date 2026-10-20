# Step 8k: register bridge-FPC tops to the RSU

The first Step 8k attempt equalized left/right offsets measured from the
*baseplate* upper y edge. The user's visual review showed why that was wrong:
the baseplates are 35 mm (A) and 32 mm (B), but the physical RSU is 19.564 mm
wide in both. Equal baseplate-relative offsets put both A FPCs 1.498 mm above
the RSU and unnecessarily moved the previously aligned B right FPC. That
attempt was not committed or accepted; its XML and validator edits were
replaced by the correction below.

The tiling builder sets each FPC upper edge to
`SiEndcapRSU_width / 2 + SiEndcapTiling[Left/Right]BridgeFPC_RSU_top_offset`.
The initial RSU-referenced revision set *both* offsets to zero. Visual review
showed that this also removed the intentional 1.28 mm left-FPC drop. The
supplied detail drawing gives 7.50 mm and 6.22 mm: their difference is
1.28 mm. The corrected offsets are now **left -1.28 mm, right 0 mm**. Thus
the right upper edge is +9.782 mm and the left upper edge is +8.502 mm in
module coordinates, for both A and B. The legacy baseplate-relative offsets
are restored unchanged and belong only to `EIC_LAS_6RSU_CORR`. The AncASIC
uses the same left-FPC y geometry. The Step 8j electronics-y correction keeps
these features at the corresponding global upper edge for negative-y modules.

The nominal FPC lengths stay 24.50 mm left and 25.78 mm right. The existing
epoxy keep-out clips only the lower FPC edges as needed: A left
[-15.000, +8.502] mm, A right [-15.000, +9.782] mm, B left
[-13.500, +8.502] mm, and B right [-13.500, +9.782] mm. The 0.5 mm x gaps
from LEC/REC, z contact with the
baseplate, sensor tiles, baseplate, module transforms, and source CSVs are
unchanged.

The Step 8d validator now asserts that the right FPC upper edge equals the
RSU upper edge and the left sits 1.28 mm below it for A/B, while retaining
its epoxy, x-gap, and z-stack checks. The
Step 8j negative-global-Y validator remains part of validation.

The previous equal-upper-edge version was visually reviewed and rejected by
the user. For the corrected 1.28 mm offset, the Step 8d and Step 8j static
validators, both shell syntax checks, and `git diff --check` passed. The
EIC-shell build/install and geometry export completed in the
`20261020T150520Z` Step 8k run. Its ROOT check found zero illegal overlaps or
extrusions at 0.01 mm, and the user accepted the visual result. The matching
ROOT report is `check_svt_fpc_top_alignment_root_step8k_20261020T150520Z.txt`.
The separate, expensive Geant4 overlap check remains unrun.

If needed, rerun `check_svt_fpc_top_alignment_visual_root_step8k.sh` in the
stable EIC shell to rebuild/export and repeat ROOT/visual checks. The next
independent check is `check_svt_fpc_top_alignment_geant_step8k.sh` for Geant4.
