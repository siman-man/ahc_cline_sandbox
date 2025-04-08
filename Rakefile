require 'open3'
require 'rake/clean'
require 'parallel'

PROBLEM_NAME = 'solver'
SEED = 0
CLEAN.include %w(data/* *.gcda *.gcov *.gcno *.png)

desc 'setup match data'
task :setup_match_data do
  Parallel.map([*1..10], in_processes: 4) do |x|
    sh("ruby simulator.rb #{x}")
  end
end


desc 'c++ file compile'
task :compile do
  sh("g++ -std=c++20 -I /Users/siman/Programming/c++/ac-library -D__NO_INLINE__ -W -Wall -Wno-sign-compare -O3 -o #{PROBLEM_NAME} #{PROBLEM_NAME}.cpp")
  # sh("g++ -std=c++20 -I /Users/siman/Programming/c++/ac-library -W -Wall -g -fsanitize=address -fno-omit-frame-pointer -Wno-sign-compare -O2 -o #{PROBLEM_NAME} #{PROBLEM_NAME}.cpp")
end

desc 'check single'
task one: [:compile] do
  sh("./#{PROBLEM_NAME} < in/%04d.txt > out/%04d.txt" % [SEED, SEED])
  sh("target/release/vis in/%04d.txt out/%04d.txt" % [SEED, SEED])
end

desc "example test"
task sample: [:compile] do
  run_test(0..9)
end

desc "production test"
task test: [:compile] do
  run_test(0..99)
end

desc "system test"
task final: [:compile] do
  run_test(1000..1999)
end

desc "system test production"
task production: [:compile] do
  run_test(0..999)
end

def run_test(seeds, options = "")
  results = Parallel.map(seeds, in_processes: 4) do |seed|
    print "\rseed = #{seed}"
    data = Open3.capture3("./#{PROBLEM_NAME} < in/%04d.txt > out/%04d.txt" % [seed, seed])
    data += Open3.capture3("target/release/vis in/%04d.txt out/%04d.txt" % [seed, seed])
    [seed, data]
  end.to_h

  File.open("result.txt", "w") do |file|
    seeds.each do |seed|
      file.puts("----- !BEGIN! ------")
      file.puts("Seed = #{seed}")

      data = results[seed]
      file.puts(data.select { |d| d.is_a?(String) }.flat_map { |d| d.split("\n") })
      file.puts("----- !END! ------")
    end
  end

  ruby "scripts/analyze.rb #{seeds.size}"
end

task default: :compile

