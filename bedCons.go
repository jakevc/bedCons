package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"os"
	"strconv"
	"strings"

	"github.com/montanaflynn/stats"

	"github.com/pbenner/gonetics"
	"github.com/pkg/errors"
)

func main() {

	bed := flag.String("bedfile", "", "BED6 format file with regions to summarize wig signal over")
	bw := flag.String("bigwig", "", "big wig file")
	mc := flag.Float64("minscore", 0.5, "the output column percentgt calcualtes the percentage of bigWig bins greater than or equal to this score")
	ofn := flag.String("outfile", "", "output filename")
	bs := flag.Int("binsize", 50, "bin size to summarize bigwig values over")
	header := flag.Bool("header", true, "print the file header")
	flag.Parse()

	if *bw == "" || *bed == "" || *ofn == "" {
		flag.Usage()
		os.Exit(1)
	}

	// query details
	binsize := *bs
	binOlap := 0
	cons := *mc
	// the bin stat
	stat := gonetics.BinSummaryStatisticsFromString("mean")

	// read the bedfile as granges
	g, err := readBed6(*bed)
	if err != nil {
		log.Fatal(err)
	}

	f, err := gonetics.OpenBigWigFile(*bw)
	defer f.Close()
	if err != nil {
		fmt.Println(errors.Wrap(err, "failed to open bigwig"))
	}

	// read the bigwig file
	reader, err := gonetics.NewBigWigReader(f)
	if err != nil {
		fmt.Println(errors.Wrap(err, "bigwig read failed"))
	}
	gr := addBinStats(g, reader, stat, binsize, binOlap, cons)
	writeOutput(gr, *ofn, *header)
}

// Import GRanges from a Bed file with 6 columns.
func readBed6(bedfile string) (gonetics.GRanges, error) {
	// open file
	f, err := os.Open(bedfile)
	if err != nil {
		return gonetics.GRanges{}, err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	// it seems that buffering the data does not increase
	// performance
	seqnames := []string{}
	from := []int{}
	to := []int{}
	name := []string{}
	score := []int{}
	strand := []byte{}

	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		if len(fields) == 0 {
			continue
		}
		// drop any track lines
		if fields[0] == "track" {
			continue
		}
		if len(fields) < 6 {
			return gonetics.GRanges{}, fmt.Errorf("ReadBed6(): Bed file must have at least 6 columns")
		}
		t1, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			return gonetics.GRanges{}, err
		}
		t2, err := strconv.ParseInt(fields[2], 10, 64)
		if err != nil {
			return gonetics.GRanges{}, err
		}
		t3, err := strconv.ParseInt(fields[4], 10, 64)
		if err != nil {
			t3 = 0
		}
		seqnames = append(seqnames, fields[0])
		from = append(from, int(t1))
		to = append(to, int(t2))
		name = append(name, fields[3])
		score = append(score, int(t3))
		if fields[5][0] == '.' {
			strand = append(strand, '*')
		} else {
			strand = append(strand, fields[5][0])
		}
	}
	g := gonetics.NewGRanges(seqnames, from, to, strand)
	g.AddMeta("name", name)
	g.AddMeta("score", score)
	return g, nil
}

// calculate the percent of an slice greater than value
func percGthan(arr []float64, value float64) float64 {
	gt := float64(0)
	n := float64(len(arr))
	for i := 0; i < len(arr); i++ {
		if arr[i] >= value {
			gt++
		}
	}
	return gt / n
}

// calculates median conservation
func addBinStats(g gonetics.GRanges, reader *gonetics.BigWigReader, binstat gonetics.BinSummaryStatistics, binsize int, binolap int, minCons float64) gonetics.GRanges {
	n := g.Length()
	medianMeta := make([]float64, n)
	meanMeta := make([]float64, n)
	percMeta := make([]float64, n)
	for i := 0; i < n; i++ {
		from := g.Ranges[i].From
		to := g.Ranges[i].To
		name := g.Seqnames[i]
		if r, _, err := reader.QuerySlice(name, from, to, binstat, binsize, binolap, math.NaN()); err != nil {
			fmt.Println(err)
		} else {
			med, err := stats.Median(r)
			if err != nil {
				fmt.Println(err)
			}

			mean, err := stats.Mean(r)
			if err != nil {
				fmt.Println(err)
			}
			medianMeta[i] = med
			meanMeta[i] = mean
			if err != nil {
				fmt.Println(err)
			}
			percMeta[i] = percGthan(r, minCons)
		}
	}
	g.AddMeta("mean", meanMeta)
	g.AddMeta("median", medianMeta)
	g.AddMeta("percentgt", percMeta)
	return g
}

func writeOutput(g gonetics.GRanges, outfile string, header bool) {
	var b bytes.Buffer

	name := g.GetMetaStr("name")
	score := g.GetMetaInt("score")
	str := "."
	mea := g.GetMetaFloat("mean")
	med := g.GetMetaFloat("median")
	per := g.GetMetaFloat("percentgt")

	if header {
		fmt.Fprintf(&b, "chrom\tstart\tend\tname\tscore\tstrand\tmean\tmedian\tpercentgt\n")
	}

	for i := 0; i < g.Length(); i++ {
		if _, err := fmt.Fprintf(&b, "%s", g.Seqnames[i]); err != nil {
			log.Fatal(err)
		}
		if _, err := fmt.Fprintf(&b, "\t%d", g.Ranges[i].From); err != nil {
			log.Fatal(err)
		}
		if _, err := fmt.Fprintf(&b, "\t%d", g.Ranges[i].To); err != nil {
			log.Fatal(err)
		}
		if _, err := fmt.Fprintf(&b, "\t%s", name[i]); err != nil {
			log.Fatal(err)
		}
		if _, err := fmt.Fprintf(&b, "\t%d", score[i]); err != nil {
			log.Fatal(err)
		}
		if _, err := fmt.Fprintf(&b, "\t%s", str); err != nil {
			log.Fatal(err)
		}
		if len(mea) > 0 {
			if _, err := fmt.Fprintf(&b, "\t%f", mea[i]); err != nil {
				log.Fatal(err)
			}
		} else {
			if _, err := fmt.Fprintf(&b, "\t%f", 0); err != nil {
				log.Fatal(err)
			}
		}
		if _, err := fmt.Fprintf(&b, "\t%f", med[i]); err != nil {
			log.Fatal(err)
		}
		if _, err := fmt.Fprintf(&b, "\t%f", per[i]); err != nil {
			log.Fatal(err)
		}
		fmt.Fprintf(&b, "\n")
	}
	ioutil.WriteFile(outfile, b.Bytes(), 0666)
}
