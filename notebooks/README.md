# notebooks/

ChCODH2 Ni-insertase in silico 스크리닝 작업용 셀 모음.

| 파일 | 용도 |
|---|---|
| `antigravity_cells_ChCODH2_screen.md` | **복사·붙여넣기용 원본.** CELL 00~33, 각 셀 상단에 번호+헤더 주석 |
| `ChCODH2_screen_antigravity.ipynb` | 위 파일에서 자동 생성한 노트북 (그대로 열어 실행해도 됨) |

## 원본 노트북(`RF2PPI_ChCODH2_screen.ipynb`)에서 바뀐 점

1. `BASE`: `/data/chcodh2_ppi` → `/mnt/af2results/mingyu/chcodh2_ppi` (서버에 `/data` 없음)
2. BL21 proteome ID: `UP000002032` → **`UP000503272`** (구조 DB/Folddisco 인덱스와 통일), Y19 `UP000034085` 추가
3. 프로테옴 서열·4-category 분류: UniProt 재다운로드 대신 **기존 GenBank 자산 재사용**, ID 체계를 GenBank로 통일
4. Folddisco 결과(`AF-<UniProt>` / `cf_<UniParc>`)를 **구조 서열 완전일치 crosswalk**로 GenBank ID에 연결해 Part 8 통합

## 그 밖의 추가

- `%%bash` 매직 제거 → 전부 파이썬 셀 (`sh` / `sh_bg` / `bg_tail` 헬퍼로 백그라운드 실행·모니터)
- segment bait 생성 (설계문서 §2 대응 b/c — apo 중간체 결합 가능성)
- gate 깊이에 따라 replicate 3회/5회 자동 결정
- full-length + segment bait 결과를 prey 단위 합집합으로 랭킹
