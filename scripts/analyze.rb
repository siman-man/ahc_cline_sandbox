require 'fileutils'
require 'pathname'

class Analyze
  TSV_LABEL = ['seed', 'score', 'time', 'd', 'n', 'max_e'].join("\t")

  def initialize
    data_list = parse
    time_stamp = Time.now.strftime("%m%d-%H%m%S")

    result_file_path = Pathname.new(File.expand_path("data/record-#{time_stamp}.tsv", Dir.pwd))
    Dir.mkdir(result_file_path.dirname) unless File.exist?(result_file_path.dirname)
    @record = File.open(result_file_path, 'w')
    result = []

    begin
      @record.puts(TSV_LABEL)

      perfect_cnt = 0
      middle_perfect_cnt = 0
      score_coutner = Hash.new(0)
      sum_score = 0.0
      seeds = []
      scores = []

      data_list.each do |data|
        data = clean_data(data)
        result << data

        seeds << data[:seed]
        scores << { seed: data[:seed], score: data[:score] }
        if data[:score] == -1.0
          sum_score += Float::INFINITY
        else
          score = data[:score]
          perfect_cnt += 1 if data[:score] == 1.0
          middle_perfect_cnt += 1 if data[:score] == 0.5
          score_coutner[(100 * score).round] += 1

          sum_score += data[:score]
        end
        @record.puts(data.values.join("\t"))
      end

      result.each do |data|
        # puts [data[:l], data[:n], data[:s]].join(',')
      end

      scores.sort_by { |h| h[:seed] }.each do |data|
        # puts Math.log(data[:score] + 1) if !ENV["MM_TUNING_FLG"]
        # puts "Seed #{data[:seed]}: #{data[:score]}" if !ENV["MM_TUNING_FLG"]
        puts data[:score] if !ENV["MM_TUNING_FLG"]
      end

      puts "%.4f" % Rational(sum_score, data_list.size)
      exec("echo %.4f > cur_score.txt" % [Rational(sum_score, data_list.size)])

      if ENV['DEBUG']
        puts "Perfect rate: %f, Middle perfect rate: %f" % [perfect_cnt / data_list.size.to_f, middle_perfect_cnt / data_list.size.to_f]

        10.times do |i|
          cnt = 0
          from = i * 10 + 1
          to = (i + 1) * 10

          from.upto(to) do |v|
            cnt += score_coutner[v]
          end

          puts "%d" % cnt
        end
      end
    ensure
      @record&.close
    end
  end 

  def clean_data(data)
    seed = data['seed']&.to_i
    score = [0, data['score']&.to_i].max
    if score == 0
      score = -1 
    end
    /(?<minute>\d+)m(?<second>(\d|\.)+)s/ =~ data['user']
    time = minute.to_f * 60 + second.to_f
    d = data['d']&.to_i
    n = data['n']&.to_i
    max_e = data['ave_e']&.to_i

    { seed: seed, score: score, time: time, d: d, n: n, max_e: max_e }
  end

  def parse
    filepath = File.expand_path('result.txt', Dir.pwd)
    data_list = []
    data = {}

    File.readlines(filepath).each do |line|
      line = clean_line(line)

      if line =~ /begin/
        data = {}
      elsif line =~ /!end!/
        data_list << data.dup
      else
        if validate(line)
          h = line2hash(line)
          data.merge!(h)
        end
      end
    end

    data_list
  end

  def line2hash(line)
    Hash[*line.split]
  end

  def clean_line(line)
    line.chomp.downcase.delete('=')
  end

  def validate(line)
    return false if line =~ /^div/
    line =~ /^(score|seed|user|time|d|n|max_e|ave_e)/
  end
end

Analyze.new
