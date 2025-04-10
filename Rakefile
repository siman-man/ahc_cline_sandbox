require 'open3'
require 'timeout'
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

# タイムアウト付きでコマンドを実行して結果を返す
def capture_with_timeout(cmd, timeout_sec)
  Open3.popen3(cmd) do |stdin, stdout, stderr, thread|
    begin
      # timeout_sec 秒以内に終了しなければTimeout::Errorになる
      Timeout.timeout(timeout_sec) do
        thread.value  # プロセス終了待ち
      end
    rescue Timeout::Error
      # 外部プロセスを kill しておかないと生き続ける可能性がある
      Process.kill('TERM', thread.pid) rescue nil
      # 必要に応じてさらに強制的に kill したい場合は KILL シグナルも
      # Process.kill('KILL', thread.pid) rescue nil
      return [
        "",
        "ERROR: Command '#{cmd}' timed out after #{timeout_sec} seconds."
      ]
    end

    # 正常終了した場合の標準出力・標準エラーを読み込む
    out = stdout.read
    err = stderr.read
    [out, err]
  end
end

def run_test(seeds, options = "")
  results = Parallel.map(seeds, in_processes: 4) do |seed|
    print "\rseed = #{seed}"

    cmd1 = "./#{PROBLEM_NAME} < in/%04d.txt > out/%04d.txt" % [seed, seed]
    cmd2 = "target/release/vis in/%04d.txt out/%04d.txt" % [seed, seed]

    out1, err1 = capture_with_timeout(cmd1, 2)
    out2, err2 = capture_with_timeout(cmd2, 2)
    [seed, [out1, err1, out2, err2]]
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

